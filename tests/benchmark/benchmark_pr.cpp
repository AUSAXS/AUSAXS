// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

// Histogram- and intensity-calculation benchmarks over the structures shipped in tests/files/.
//
// Run (from the repository root, so the tests/files/ paths resolve):
//     cmake --build build --target benchmarks -j8
//     ./build/tests/benchmark/bin/benchmark_pr "[benchmark]~[large]" -r xml::out=results.xml
//     python3 tests/benchmark/plot_benchmark.py results.xml
//
// The [large] tag guards A2M_native (43k atoms): a single Grid histogram there is minutes, so a
// default-sample Catch2 run of it is hours. Give it its own run with --benchmark-samples 3.

#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <data/Body.h>
#include <data/Molecule.h>
#include <data/symmetry/PointSymmetry.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFExplicit.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFGrid.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogram.h>
#include <io/File.h>
#include <math/Vector3.h>
#include <rigidbody/BodySplitter.h>
#include <settings/All.h>
#include <utility/Random.h>

#include <array>

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

// Every structure here ships in tests/files/, so the benchmarks run from a clean clone.
struct MolSpec { const char* pdb; const char* label; };
static constexpr std::array<MolSpec, 5> bench_molecules = {{
    {.pdb="tests/files/2epe.pdb",        .label="2epe"},
    {.pdb="tests/files/LAR1-2.pdb",      .label="LAR1-2"},
    {.pdb="tests/files/SASDJQ4.pdb",     .label="SASDJQ4"},
    {.pdb="tests/files/SASDJG5.pdb",     .label="SASDJG5"},
    {.pdb="tests/files/168l.pdb",        .label="168l"},
}};
static constexpr MolSpec large_molecule = {.pdb="tests/files/A2M_native.pdb", .label="A2M_native"};

// Benchmarks measure the shipped default configuration: implicit hydrogens stay on, since that is
// what `ausaxs fit` runs and it is what decides how many form-factor slots are active.
static void bench_settings() {
    settings::general::verbose = false;
    settings::general::warnings = false;
    settings::molecule::allow_unknown_residues = true;
    random::set_seed(0); // hydration is otherwise unseeded, and water count feeds every histogram below
}

// Section label carrying the atom count, which plot_benchmark.py reads back off the XML.
static std::string section_label(const MolSpec& spec) {
    data::Molecule mol(spec.pdb);
    return std::string(spec.label) + " (" + std::to_string(mol.size_atom()) + " atoms)";
}

static void histogram_benchmarks(const MolSpec& spec) {
    SECTION(section_label(spec)) {
        BENCHMARK_ADVANCED("Simple") (Catch::Benchmark::Chronometer meter) {
            data::Molecule mol(spec.pdb);
            mol.set_histogram_manager(settings::hist::HistogramManagerChoice::HistogramManagerMT);
            meter.measure([&] { return mol.get_histogram(); });
        };

        BENCHMARK_ADVANCED("Fraser") (Catch::Benchmark::Chronometer meter) {
            data::Molecule mol(spec.pdb);
            mol.set_histogram_manager(settings::hist::HistogramManagerChoice::HistogramManagerMTFFExplicit);
            meter.measure([&] { return mol.get_histogram(); });
        };

        BENCHMARK_ADVANCED("Grid") (Catch::Benchmark::Chronometer meter) {
            data::Molecule mol(spec.pdb);
            mol.set_histogram_manager(settings::hist::HistogramManagerChoice::HistogramManagerMTFFGrid);
            meter.measure([&] { return mol.get_histogram(); });
        };
    }
}

static void intensity_benchmarks(const MolSpec& spec) {
    SECTION(section_label(spec)) {
        // Simple has no sinqd cache — every call IS the full Debye summation.
        BENCHMARK_ADVANCED("Simple") (Catch::Benchmark::Chronometer meter) {
            data::Molecule mol(spec.pdb);
            mol.set_histogram_manager(settings::hist::HistogramManagerChoice::HistogramManagerMT);
            auto hist = mol.get_histogram();
            meter.measure([&] { return hist->debye_transform(); });
        };

        // Fraser: sinqd cache warm — measures only intensity-profile assembly from cached sinqd values.
        BENCHMARK_ADVANCED("Fraser (sinqd warm)") (Catch::Benchmark::Chronometer meter) {
            data::Molecule mol(spec.pdb);
            mol.set_histogram_manager(settings::hist::HistogramManagerChoice::HistogramManagerMTFFExplicit);
            auto ptr = std::unique_ptr<hist::CompositeDistanceHistogramFFExplicit>(
                dynamic_cast<hist::CompositeDistanceHistogramFFExplicit*>(mol.get_histogram().release()));
            BenchFFExplicit bench(std::move(*ptr));
            bench.debye_transform(); // prime: allocate containers and fill sinqd cache
            meter.measure([&] { return bench.debye_transform(); });
        };

        // Fraser: sinqd cache cold — measures the full Debye inner-product pass + assembly.
        BENCHMARK_ADVANCED("Fraser (sinqd cold)") (Catch::Benchmark::Chronometer meter) {
            data::Molecule mol(spec.pdb);
            mol.set_histogram_manager(settings::hist::HistogramManagerChoice::HistogramManagerMTFFExplicit);
            auto ptr = std::unique_ptr<hist::CompositeDistanceHistogramFFExplicit>(
                dynamic_cast<hist::CompositeDistanceHistogramFFExplicit*>(mol.get_histogram().release()));
            BenchFFExplicit bench(std::move(*ptr));
            bench.debye_transform(); // prime: allocate containers
            meter.measure([&] {
                bench.invalidate_sinqd();
                return bench.debye_transform();
            });
        };

        // Grid: sinqd cache warm.
        BENCHMARK_ADVANCED("Grid (sinqd warm)") (Catch::Benchmark::Chronometer meter) {
            data::Molecule mol(spec.pdb);
            mol.set_histogram_manager(settings::hist::HistogramManagerChoice::HistogramManagerMTFFGrid);
            auto ptr = std::unique_ptr<hist::CompositeDistanceHistogramFFGrid>(
                dynamic_cast<hist::CompositeDistanceHistogramFFGrid*>(mol.get_histogram().release()));
            BenchFFGrid bench(std::move(*ptr));
            bench.debye_transform(); // prime
            meter.measure([&] { return bench.debye_transform(); });
        };

        // Grid: sinqd cache cold.
        BENCHMARK_ADVANCED("Grid (sinqd cold)") (Catch::Benchmark::Chronometer meter) {
            data::Molecule mol(spec.pdb);
            mol.set_histogram_manager(settings::hist::HistogramManagerChoice::HistogramManagerMTFFGrid);
            auto ptr = std::unique_ptr<hist::CompositeDistanceHistogramFFGrid>(
                dynamic_cast<hist::CompositeDistanceHistogramFFGrid*>(mol.get_histogram().release()));
            BenchFFGrid bench(std::move(*ptr));
            bench.debye_transform(); // prime
            meter.measure([&] {
                bench.invalidate_sinqd();
                return bench.debye_transform();
            });
        };
    }
}

TEST_CASE("Distance calculation benchmark: real molecules", "[.][benchmark]") {
    bench_settings();
    histogram_benchmarks(bench_molecules[GENERATE(0, 1, 2, 3, 4)]);
}

TEST_CASE("Distance calculation benchmark: large molecule", "[.][benchmark][large]") {
    bench_settings();
    histogram_benchmarks(large_molecule);
}

TEST_CASE("Intensity calculation benchmark: debye_transform only", "[.][benchmark]") {
    bench_settings();
    intensity_benchmarks(bench_molecules[GENERATE(0, 1, 2, 3, 4)]);
}

TEST_CASE("Intensity calculation benchmark: large molecule", "[.][benchmark][large]") {
    bench_settings();
    intensity_benchmarks(large_molecule);
}

TEST_CASE("Symmetry histogram benchmark", "[.][benchmark]") {
    bench_settings();

    SECTION("SASDJG5 (dimer symmetry)") {
        auto mol = rigidbody::BodySplitter::split(io::File("tests/files/SASDJG5_single.pdb"));
        auto mol_full = rigidbody::BodySplitter::split(io::File("tests/files/SASDJG5.pdb"));
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
                auto* s = mol.get_body(0).symmetry().get(0);
                static_cast<symmetry::PointSymmetry*>(s)->translation += Vector3<double>{0.1, 0, 0};
                return mol.get_histogram()->debye_transform();
            });
        };

        BENCHMARK_ADVANCED("PartialHistogramSymmetryManagerMT") (Catch::Benchmark::Chronometer meter) {
            mol.set_histogram_manager(settings::hist::HistogramManagerChoice::PartialHistogramSymmetryManagerMT);
            [[maybe_unused]] auto warmup = mol.get_histogram();
            meter.measure([&] {
                auto* s = mol.get_body(0).symmetry().get(0);
                static_cast<symmetry::PointSymmetry*>(s)->translation += Vector3<double>{0.1, 0, 0};
                return mol.get_histogram()->debye_transform();
            });
        };
    }
}

TEST_CASE("Partial histogram benchmark: rigidbody body translation", "[.][benchmark]") {
    bench_settings();

    struct SplitSpec { const char* pdb; const char* label; std::vector<int> splits; };
    auto spec = GENERATE(
        SplitSpec{"tests/files/SASDJG5.pdb", "SASDJG5 (chain split, 2 bodies)",  {}},
        SplitSpec{"tests/files/LAR1-2.pdb",  "LAR1-2 (index split, 4 bodies)",   {50, 100, 150}},
        SplitSpec{"tests/files/168l.pdb",    "168l (index split, 3 bodies)",     {50, 100}}
    );

    SECTION(spec.label) {
        auto make = [&] {
            return spec.splits.empty()
                ? rigidbody::BodySplitter::split(io::File(spec.pdb))
                : rigidbody::BodySplitter::split(io::File(spec.pdb), spec.splits);
        };
        auto mol = make();
        INFO(mol.size_atom() << " atoms, " << mol.size_body() << " bodies");

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
