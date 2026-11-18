// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

// Grid cell width benchmarks: what the default hydration grid width costs on small structures.
//
// Run (from the repository root, so the tests/files/ paths resolve):
//     cmake --build build --target benchmarks
//     ./build/tests/benchmark/bin/benchmark_grid_width "[benchmark]" --benchmark-samples 20
//
// Each stage is timed at 1.0, 1.5 and 1.75 A:
//   grid          - allocating the grid and marking the atom centres (Molecule::create_grid)
//   grid + volume - the same, plus expanding every atom's volume into it, which is the part that scales as width^-3
//   hydration     - one generate_new_hydration() on an existing grid
//   end-to-end    - what `ausaxs fit` does for a structure: load, grid, hydrate, histogram, fit

#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_range.hpp>

#include <data/Molecule.h>
#include <dataset/SimpleDataset.h>
#include <fitter/SmartFitter.h>
#include <form_factor/lookup/FormFactorManager.h>
#include <grid/Grid.h>
#include <hist/scattering_model/IScatteringModel.h> // NOLINT
#include <io/File.h>
#include <settings/All.h>

#include <array>
#include <string>

using namespace ausaxs;

namespace {
    struct Spec {const char* pdb; const char* dat; const char* label;};
    constexpr std::array<Spec, 4> structures = {{
        {.pdb="tests/files/2epe.pdb",    .dat="tests/files/2epe.dat",    .label="2epe"},
        {.pdb="tests/files/LAR1-2.pdb",  .dat="tests/files/LAR1-2.dat",  .label="LAR1-2"},
        {.pdb="tests/files/SASDJQ4.pdb", .dat="tests/files/SASDJQ4.dat", .label="SASDJQ4"},
        {.pdb="tests/files/SASDJG5.pdb", .dat="tests/files/SASDJG5.dat", .label="SASDJG5"},
    }};

    void bench_settings(double width) {
        settings::general::verbose = false;
        settings::general::warnings = false;
        settings::grid::cell_width = width;
    }
}

TEST_CASE("Grid width benchmark", "[.][benchmark]") {
    double width = GENERATE(1.0, 1.5, 1.75);
    auto spec = GENERATE(from_range(structures));
    bench_settings(width);

    SECTION(std::string(spec.label) + " @ " + std::to_string(width).substr(0, 4) + " A") {
        BENCHMARK_ADVANCED("grid") (Catch::Benchmark::Chronometer meter) {
            data::Molecule mol(spec.pdb);
            meter.measure([&] {mol.clear_grid(); return mol.create_grid();});
        };

        BENCHMARK_ADVANCED("grid + volume") (Catch::Benchmark::Chronometer meter) {
            data::Molecule mol(spec.pdb);
            meter.measure([&] {mol.clear_grid(); auto* grid = mol.create_grid(); grid->expand_volume(); return grid->get_volume_bins();});
        };

        BENCHMARK_ADVANCED("hydration") (Catch::Benchmark::Chronometer meter) {
            data::Molecule mol(spec.pdb);
            mol.generate_new_hydration();
            meter.measure([&] {mol.generate_new_hydration(); return mol.size_water();});
        };

        BENCHMARK("end-to-end") {
            data::Molecule mol(spec.pdb);
            mol.generate_new_hydration();
            form_factor::manager::use_form_factors(mol);
            SimpleDataset data(spec.dat);
            fitter::SmartFitter fitter(data, mol.get_histogram());
            return fitter.fit()->fval;
        };
    }
}
