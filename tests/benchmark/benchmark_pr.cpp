#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>

#include <hist/histogram_manager/HistogramManagerFactory.h>
#include <hist/histogram_manager/IPartialHistogramManager.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogram.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFExplicit.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFGrid.h>
#include <data/symmetry/PointSymmetry.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <io/File.h>
#include <math/Vector3.h>
#include <rigidbody/BodySplitter.h>
#include <settings/All.h>

using namespace ausaxs;

// Expose cache.sinqd.valid for Fraser and Grid intensity calculators.
struct BenchFFExplicit : public hist::CompositeDistanceHistogramFFExplicit {
    BenchFFExplicit(hist::CompositeDistanceHistogramFFExplicit&& o)
        : hist::CompositeDistanceHistogramFFExplicit(std::move(o)) {}
    void invalidate_sinqd() { cache.sinqd.valid = false; }
};

struct BenchFFGrid : public hist::CompositeDistanceHistogramFFGrid {
    BenchFFGrid(hist::CompositeDistanceHistogramFFGrid&& o)
        : hist::CompositeDistanceHistogramFFGrid(std::move(o)) {}
    void invalidate_sinqd() { cache.sinqd.valid = false; }
};

struct MolSpec { const char* pdb; const char* label; };
static constexpr MolSpec bench_molecules[] = {
    {"data/SASDE35/SASDE35.pdb",       "SASDE35"},
    {"data/SASDJG5/SASDJG5.pdb",       "SASDJG5"},
    {"data/SASDME4/SASDME4.pdb",       "SASDME4"},
    {"data/SASDDD3/SASDDD3.pdb",       "SASDDD3"},
    {"data/SASDA45/SASDA45.pdb",       "SASDA45"},
    {"data/A2M_native/A2M_native.pdb", "A2M_native"},
};

// Stripped output paths — same directory, suffix "_bench_stripped.pdb".
static std::string stripped_path(const MolSpec& s) {
    std::string p = s.pdb;
    auto slash = p.rfind('/');
    auto dir   = p.substr(0, slash + 1);
    return dir + s.label + "_bench_stripped.pdb";
}

// Atom counts filled in by the setup test case and consumed by test 1.
static std::vector<std::size_t> g_atom_counts;

TEST_CASE("Benchmark setup: write stripped PDB files", "[.][benchmark]") {
    settings::general::verbose = false;
    settings::molecule::allow_unknown_residues = true;
    settings::molecule::implicit_hydrogens = false;
    g_atom_counts.clear();
    for (const auto& s : bench_molecules) {
        data::Molecule mol(s.pdb);
        g_atom_counts.push_back(mol.size_atom());
        mol.save(io::File(stripped_path(s)));
    }
    SUCCEED(); // no assertions — just setup
}

TEST_CASE("Distance calculation benchmark: real molecules", "[.][benchmark]") {
    settings::general::verbose = false;
    settings::molecule::allow_unknown_residues = true;
    settings::molecule::implicit_hydrogens = false;

    auto idx = GENERATE(0, 1, 2, 3, 4, 5);
    const MolSpec& spec = bench_molecules[idx];
    std::string path  = stripped_path(spec);
    std::size_t n_atoms = g_atom_counts.empty() ? 0 : g_atom_counts[idx];
    std::string section_label = std::string(spec.label)
        + (n_atoms ? " (" + std::to_string(n_atoms) + " atoms)" : "");

    SECTION(section_label) {
        BENCHMARK_ADVANCED("Simple") (Catch::Benchmark::Chronometer meter) {
            meter.measure([&] {
                data::Molecule mol(path);
                mol.set_histogram_manager(settings::hist::HistogramManagerChoice::HistogramManagerMT);
                return mol.get_histogram()->debye_transform();
            });
        };

        BENCHMARK_ADVANCED("Fraser") (Catch::Benchmark::Chronometer meter) {
            meter.measure([&] {
                data::Molecule mol(path);
                mol.set_histogram_manager(settings::hist::HistogramManagerChoice::HistogramManagerMTFFExplicit);
                return mol.get_histogram()->debye_transform();
            });
        };

        BENCHMARK_ADVANCED("Grid") (Catch::Benchmark::Chronometer meter) {
            meter.measure([&] {
                data::Molecule mol(path);
                mol.set_histogram_manager(settings::hist::HistogramManagerChoice::HistogramManagerMTFFGrid);
                return mol.get_histogram()->debye_transform();
            });
        };
    }
}

TEST_CASE("Intensity calculation benchmark: debye_transform only", "[.][benchmark]") {
    settings::general::verbose = false;
    settings::molecule::allow_unknown_residues = true;
    settings::molecule::implicit_hydrogens = false;

    auto idx = GENERATE(0, 1, 2, 3, 4, 5);
    const MolSpec& spec = bench_molecules[idx];
    std::string path  = stripped_path(spec);
    std::size_t n_atoms = g_atom_counts.empty() ? 0 : g_atom_counts[idx];
    std::string section_label = std::string(spec.label)
        + (n_atoms ? " (" + std::to_string(n_atoms) + " atoms)" : "");

    SECTION(section_label) {
        // Simple has no sinqd cache — every call IS the full Debye summation.
        BENCHMARK_ADVANCED("Simple") (Catch::Benchmark::Chronometer meter) {
            data::Molecule mol(path);
            mol.set_histogram_manager(settings::hist::HistogramManagerChoice::HistogramManagerMT);
            auto hist = mol.get_histogram();
            meter.measure([&] {
                return hist->debye_transform();
            });
        };

        // Fraser: sinqd cache warm — measures only intensity-profile assembly from cached sinqd values.
        BENCHMARK_ADVANCED("Fraser (sinqd warm)") (Catch::Benchmark::Chronometer meter) {
            data::Molecule mol(path);
            mol.set_histogram_manager(settings::hist::HistogramManagerChoice::HistogramManagerMTFFExplicit);
            auto ptr = std::unique_ptr<hist::CompositeDistanceHistogramFFExplicit>(
                dynamic_cast<hist::CompositeDistanceHistogramFFExplicit*>(mol.get_histogram().release()));
            BenchFFExplicit bench(std::move(*ptr));
            bench.debye_transform(); // prime: allocate containers and fill sinqd cache
            meter.measure([&] {
                return bench.debye_transform();
            });
        };

        // Fraser: sinqd cache cold — measures the full Debye inner-product pass + assembly.
        BENCHMARK_ADVANCED("Fraser (sinqd cold)") (Catch::Benchmark::Chronometer meter) {
            data::Molecule mol(path);
            mol.set_histogram_manager(settings::hist::HistogramManagerChoice::HistogramManagerMTFFExplicit);
            auto ptr = std::unique_ptr<hist::CompositeDistanceHistogramFFExplicit>(
                dynamic_cast<hist::CompositeDistanceHistogramFFExplicit*>(mol.get_histogram().release()));
            BenchFFExplicit bench(std::move(*ptr));
            bench.debye_transform(); // prime: allocate containers
            bench.invalidate_sinqd();
            meter.measure([&] {
                bench.invalidate_sinqd();
                return bench.debye_transform();
            });
        };

        // Grid: sinqd cache warm.
        BENCHMARK_ADVANCED("Grid (sinqd warm)") (Catch::Benchmark::Chronometer meter) {
            data::Molecule mol(path);
            mol.set_histogram_manager(settings::hist::HistogramManagerChoice::HistogramManagerMTFFGrid);
            auto ptr = std::unique_ptr<hist::CompositeDistanceHistogramFFGrid>(
                dynamic_cast<hist::CompositeDistanceHistogramFFGrid*>(mol.get_histogram().release()));
            BenchFFGrid bench(std::move(*ptr));
            bench.debye_transform(); // prime
            meter.measure([&] {
                return bench.debye_transform();
            });
        };

        // Grid: sinqd cache cold.
        BENCHMARK_ADVANCED("Grid (sinqd cold)") (Catch::Benchmark::Chronometer meter) {
            data::Molecule mol(path);
            mol.set_histogram_manager(settings::hist::HistogramManagerChoice::HistogramManagerMTFFGrid);
            auto ptr = std::unique_ptr<hist::CompositeDistanceHistogramFFGrid>(
                dynamic_cast<hist::CompositeDistanceHistogramFFGrid*>(mol.get_histogram().release()));
            BenchFFGrid bench(std::move(*ptr));
            bench.debye_transform(); // prime
            bench.invalidate_sinqd();
            meter.measure([&] {
                bench.invalidate_sinqd();
                return bench.debye_transform();
            });
        };
    }
}

TEST_CASE("Symmetry histogram benchmark", "[.][benchmark]") {
    settings::general::verbose = false;
    settings::molecule::allow_unknown_residues = true;
    settings::molecule::implicit_hydrogens = false;

    SECTION("SASDJG5 (dimer symmetry)") {
        auto mol = rigidbody::BodySplitter::split("data/SASDJG5/SASDJG5_single.pdb");
        auto mol_full = rigidbody::BodySplitter::split("data/SASDJG5/SASDJG5.pdb");
        mol.get_body(0).symmetry().add(symmetry::type::p2);
        INFO(mol.size_atom() << " atoms");

        BENCHMARK_ADVANCED("HistogramManagerMT") (Catch::Benchmark::Chronometer meter) {
            mol_full.set_histogram_manager(settings::hist::HistogramManagerChoice::HistogramManagerMT);
            meter.measure([&] {
                mol_full.get_body(0).translate(Vector3<double>(0.1, 0, 0));
                return mol_full.get_histogram()->debye_transform();
            });
        };

        BENCHMARK_ADVANCED("PartialHistogramManagerMT") (Catch::Benchmark::Chronometer meter) {
            mol_full.set_histogram_manager(settings::hist::HistogramManagerChoice::PartialHistogramManagerMT);
            [[maybe_unused]] auto warmup = mol_full.get_histogram();
            meter.measure([&] {
                mol_full.get_body(0).translate(Vector3<double>(0.1, 0, 0));
                return mol_full.get_histogram()->debye_transform();
            });
        };

        BENCHMARK_ADVANCED("HistogramSymmetryManagerMT") (Catch::Benchmark::Chronometer meter) {
            mol.set_histogram_manager(settings::hist::HistogramManagerChoice::HistogramSymmetryManagerMT);
            meter.measure([&] {
                auto s = mol.get_body(0).symmetry().get(0);
                static_cast<symmetry::PointSymmetry*>(s)->translation += Vector3<double>{0.1, 0, 0};
                return mol.get_histogram()->debye_transform();
            });
        };

        BENCHMARK_ADVANCED("PartialHistogramSymmetryManagerMT") (Catch::Benchmark::Chronometer meter) {
            mol.set_histogram_manager(settings::hist::HistogramManagerChoice::PartialHistogramSymmetryManagerMT);
            [[maybe_unused]] auto warmup = mol.get_histogram();
            meter.measure([&] {
                auto s = mol.get_body(0).symmetry().get(0);
                static_cast<symmetry::PointSymmetry*>(s)->translation += Vector3<double>{0.1, 0, 0};
                return mol.get_histogram()->debye_transform();
            });
        };
    }    
}

TEST_CASE("Partial histogram benchmark: rigidbody body translation", "[.][benchmark]") {
    settings::general::verbose = false;
    settings::molecule::allow_unknown_residues = true;
    settings::molecule::implicit_hydrogens = false;

    SECTION("SASDJG5 (chain split, 2 bodies)") {
        auto mol = rigidbody::BodySplitter::split("data/SASDJG5/SASDJG5.pdb");
        INFO(mol.size_atom() << " atoms");

        BENCHMARK_ADVANCED("Full (baseline)") (Catch::Benchmark::Chronometer meter) {
            mol.set_histogram_manager(settings::hist::HistogramManagerChoice::HistogramManagerMT);
            meter.measure([&] {
                mol.get_body(0).translate(Vector3<double>(0.1, 0, 0));
                return mol.get_histogram()->debye_transform();
            });
        };

        BENCHMARK_ADVANCED("PartialMT") (Catch::Benchmark::Chronometer meter) {
            mol.set_histogram_manager(settings::hist::HistogramManagerChoice::PartialHistogramManagerMT);
            [[maybe_unused]] auto warmup = mol.get_histogram();
            meter.measure([&] {
                mol.get_body(0).translate(Vector3<double>(0.1, 0, 0));
                return mol.get_histogram()->debye_transform();
            });
        };
    }

    SECTION("LAR1-2-full (index split, 4 bodies)") {
        auto mol = rigidbody::BodySplitter::split("data/rigidbody/LAR1-2-full/LAR1-2-full.pdb", {98, 196, 291});
        INFO(mol.size_atom() << " atoms");

        BENCHMARK_ADVANCED("Full (baseline)") (Catch::Benchmark::Chronometer meter) {
            mol.set_histogram_manager(settings::hist::HistogramManagerChoice::HistogramManagerMT);
            meter.measure([&] {
                mol.get_body(0).translate(Vector3<double>(0.1, 0, 0));
                return mol.get_histogram()->debye_transform();
            });
        };

        BENCHMARK_ADVANCED("PartialMT") (Catch::Benchmark::Chronometer meter) {
            mol.set_histogram_manager(settings::hist::HistogramManagerChoice::PartialHistogramManagerMT);
            [[maybe_unused]] auto warmup = mol.get_histogram();
            meter.measure([&] {
                mol.get_body(0).translate(Vector3<double>(0.1, 0, 0));
                return mol.get_histogram()->debye_transform();
            });
        };
    }

    SECTION("urateox (index split, 3 bodies)") {
        auto mol = rigidbody::BodySplitter::split("data/rigidbody/urateox/urateox_stripped.pdb", {100, 200});
        INFO(mol.size_atom() << " atoms");

        BENCHMARK_ADVANCED("Full (baseline)") (Catch::Benchmark::Chronometer meter) {
            mol.set_histogram_manager(settings::hist::HistogramManagerChoice::HistogramManagerMT);
            meter.measure([&] {
                mol.get_body(0).translate(Vector3<double>(0.1, 0, 0));
                return mol.get_histogram()->debye_transform();
            });
        };

        BENCHMARK_ADVANCED("PartialMT") (Catch::Benchmark::Chronometer meter) {
            mol.set_histogram_manager(settings::hist::HistogramManagerChoice::PartialHistogramManagerMT);
            [[maybe_unused]] auto warmup = mol.get_histogram();
            meter.measure([&] {
                mol.get_body(0).translate(Vector3<double>(0.1, 0, 0));
                return mol.get_histogram()->debye_transform();
            });
        };
    }

}
