#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>

#include <hist/histogram_manager/HistogramManagerFactory.h>
#include <hist/histogram_manager/IPartialHistogramManager.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogram.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <io/File.h>
#include <math/Vector3.h>
#include <rigidbody/BodySplitter.h>
#include <settings/All.h>

using namespace ausaxs;

struct MolSpec { const char* pdb; const char* label; };

// Source PDBs for the full-recalculation benchmark (ausaxs_bench dataset).
// Each entry is loaded once during setup to write a stripped file and record
// the real atom count. The stripped files are reusable by external tools
// (CRYSOL, FoXS, Pepsi-SAXS).
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

// ────────────────────────────────────────────────────────────────────────────
// 0. Setup: load each molecule once, record atom counts, write stripped PDBs.
//    Run this before the benchmark tests by including "[benchmark]" in the filter.
// ────────────────────────────────────────────────────────────────────────────
TEST_CASE("Benchmark setup: write stripped PDB files", "[.][benchmark]") {
    settings::general::verbose = false;
    settings::molecule::allow_unknown_residues = true;
    g_atom_counts.clear();
    for (const auto& s : bench_molecules) {
        data::Molecule mol(s.pdb);
        g_atom_counts.push_back(mol.size_atom());
        mol.save(io::File(stripped_path(s)));
    }
    SUCCEED(); // no assertions — just setup
}

// ────────────────────────────────────────────────────────────────────────────
// 1. Basic histogram manager benchmark — full recalculation on every call.
//    Tests Simple (HistogramManagerMT), Fraser (FFExplicit), and Grid (FFGrid)
//    against the full range of real molecules from the ausaxs_bench dataset.
//    File loading is included in the timed section.
// ────────────────────────────────────────────────────────────────────────────
TEST_CASE("Distance calculation benchmark: real molecules", "[.][benchmark]") {
    settings::general::verbose = false;
    settings::molecule::allow_unknown_residues = true;

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

// ────────────────────────────────────────────────────────────────────────────
// 2. Partial histogram manager benchmark — incremental update cost.
//    Models the inner loop of rigidbody optimisation: one body is translated
//    and only the affected partial histograms are recomputed.
//
//    Each BENCHMARK_ADVANCED re-initialises the manager and does one full
//    warm-up calculation before timing begins, so only the incremental update
//    cost is measured.  Two variants are compared:
//      • Full (baseline) – HistogramManagerMT, full recalculation every call
//      • PartialMT       – PartialHistogramManagerMT, multi-threaded incremental
// ────────────────────────────────────────────────────────────────────────────
TEST_CASE("Partial histogram benchmark: rigidbody body translation", "[.][benchmark]") {
    settings::general::verbose = false;
    settings::molecule::allow_unknown_residues = true;

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
