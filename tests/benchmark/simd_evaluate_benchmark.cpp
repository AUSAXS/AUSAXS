// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

/**
 * @brief Benchmarks for the distance histogram calculation.
 *
 * Measures the full mol.get_histogram() throughput using real atomic data
 * from 6lyz.pdb. Compiled in multiple SIMD variants via CMake so each binary
 * exercises exactly one dispatch level: scalar/SSE2/AVX/AVX-512.
 *
 * Run with:
 *   ./bin/benchmark_histogram_<variant> [bench]
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>

#include <data/Molecule.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogram.h>
#include <settings/All.h>

using namespace ausaxs;

TEST_CASE("histogram: distance calculation throughput", "[bench]") {
    settings::general::verbose = false;
    settings::exv::exv_method = settings::exv::ExvMethod::Simple;
    settings::hist::weighted_bins = true;
    data::Molecule mol("tests/files/6lyz.pdb");

    BENCHMARK_ADVANCED("get_histogram")(Catch::Benchmark::Chronometer meter) {
        std::unique_ptr<hist::ICompositeDistanceHistogram> hist;
        meter.measure([&] {
            hist = mol.get_histogram();
        });
        REQUIRE(hist != nullptr);
    };
}
