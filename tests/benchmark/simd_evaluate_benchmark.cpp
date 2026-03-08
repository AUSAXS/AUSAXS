// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

/**
 * @brief Microbenchmarks for the SIMD evaluate_rounded implementations on
 *        CompactCoordinatesXYZW and CompactCoordinatesXYZFF.
 *
 * Measures the throughput of computing binned distances from one reference
 * atom to an array of N atoms, using chunks of 4 (SSE quad) or 8 (SSE/AVX
 * octo), compared to the scalar baseline.
 *
 * Run with:
 *   ./bin/benchmark_simd_evaluate [bench]
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>

#include <hist/detail/data/CompactCoordinatesXYZW.h>
#include <hist/detail/data/CompactCoordinatesXYZFF.h>
#include <math/Vector3.h>
#include <form_factor/FormFactorType.h>

#include <random>
#include <vector>

using namespace ausaxs;
using namespace ausaxs::hist::detail;

// Subclass wrappers to expose the protected scalar/SSE/AVX methods
template<bool vbw>
struct XYZWBench : CompactCoordinatesXYZW<vbw> {
    using CompactCoordinatesXYZW<vbw>::CompactCoordinatesXYZW;

    xyzw::QuadEvaluatedResultRounded scalar4(
        const CompactCoordinatesXYZW<vbw>& v1, const CompactCoordinatesXYZW<vbw>& v2,
        const CompactCoordinatesXYZW<vbw>& v3, const CompactCoordinatesXYZW<vbw>& v4) const {
        return CompactCoordinatesXYZW<vbw>::evaluate_rounded_scalar(v1, v2, v3, v4);
    }
    xyzw::OctoEvaluatedResultRounded scalar8(
        const CompactCoordinatesXYZW<vbw>& v1, const CompactCoordinatesXYZW<vbw>& v2,
        const CompactCoordinatesXYZW<vbw>& v3, const CompactCoordinatesXYZW<vbw>& v4,
        const CompactCoordinatesXYZW<vbw>& v5, const CompactCoordinatesXYZW<vbw>& v6,
        const CompactCoordinatesXYZW<vbw>& v7, const CompactCoordinatesXYZW<vbw>& v8) const {
        return CompactCoordinatesXYZW<vbw>::evaluate_rounded_scalar(v1, v2, v3, v4, v5, v6, v7, v8);
    }

    #if defined __SSE2__
    xyzw::QuadEvaluatedResultRounded sse4(
        const CompactCoordinatesXYZW<vbw>& v1, const CompactCoordinatesXYZW<vbw>& v2,
        const CompactCoordinatesXYZW<vbw>& v3, const CompactCoordinatesXYZW<vbw>& v4) const {
        return CompactCoordinatesXYZW<vbw>::evaluate_rounded_sse(v1, v2, v3, v4);
    }
    xyzw::OctoEvaluatedResultRounded sse8(
        const CompactCoordinatesXYZW<vbw>& v1, const CompactCoordinatesXYZW<vbw>& v2,
        const CompactCoordinatesXYZW<vbw>& v3, const CompactCoordinatesXYZW<vbw>& v4,
        const CompactCoordinatesXYZW<vbw>& v5, const CompactCoordinatesXYZW<vbw>& v6,
        const CompactCoordinatesXYZW<vbw>& v7, const CompactCoordinatesXYZW<vbw>& v8) const {
        return CompactCoordinatesXYZW<vbw>::evaluate_rounded_sse(v1, v2, v3, v4, v5, v6, v7, v8);
    }
    #endif

    #if defined __AVX__
    xyzw::OctoEvaluatedResultRounded avx8(
        const CompactCoordinatesXYZW<vbw>& v1, const CompactCoordinatesXYZW<vbw>& v2,
        const CompactCoordinatesXYZW<vbw>& v3, const CompactCoordinatesXYZW<vbw>& v4,
        const CompactCoordinatesXYZW<vbw>& v5, const CompactCoordinatesXYZW<vbw>& v6,
        const CompactCoordinatesXYZW<vbw>& v7, const CompactCoordinatesXYZW<vbw>& v8) const {
        return CompactCoordinatesXYZW<vbw>::evaluate_rounded_avx(v1, v2, v3, v4, v5, v6, v7, v8);
    }
    #endif
};

template<bool vbw>
struct XYZFFBench : CompactCoordinatesXYZFF<vbw> {
    using CompactCoordinatesXYZFF<vbw>::CompactCoordinatesXYZFF;

    xyzff::QuadEvaluatedResultRounded scalar4(
        const CompactCoordinatesXYZFF<vbw>& v1, const CompactCoordinatesXYZFF<vbw>& v2,
        const CompactCoordinatesXYZFF<vbw>& v3, const CompactCoordinatesXYZFF<vbw>& v4) const {
        return CompactCoordinatesXYZFF<vbw>::evaluate_rounded_scalar(v1, v2, v3, v4);
    }
    xyzff::OctoEvaluatedResultRounded scalar8(
        const CompactCoordinatesXYZFF<vbw>& v1, const CompactCoordinatesXYZFF<vbw>& v2,
        const CompactCoordinatesXYZFF<vbw>& v3, const CompactCoordinatesXYZFF<vbw>& v4,
        const CompactCoordinatesXYZFF<vbw>& v5, const CompactCoordinatesXYZFF<vbw>& v6,
        const CompactCoordinatesXYZFF<vbw>& v7, const CompactCoordinatesXYZFF<vbw>& v8) const {
        return CompactCoordinatesXYZFF<vbw>::evaluate_rounded_scalar(v1, v2, v3, v4, v5, v6, v7, v8);
    }

    #if defined __SSE2__
    xyzff::QuadEvaluatedResultRounded sse4(
        const CompactCoordinatesXYZFF<vbw>& v1, const CompactCoordinatesXYZFF<vbw>& v2,
        const CompactCoordinatesXYZFF<vbw>& v3, const CompactCoordinatesXYZFF<vbw>& v4) const {
        return CompactCoordinatesXYZFF<vbw>::evaluate_rounded_sse(v1, v2, v3, v4);
    }
    xyzff::OctoEvaluatedResultRounded sse8(
        const CompactCoordinatesXYZFF<vbw>& v1, const CompactCoordinatesXYZFF<vbw>& v2,
        const CompactCoordinatesXYZFF<vbw>& v3, const CompactCoordinatesXYZFF<vbw>& v4,
        const CompactCoordinatesXYZFF<vbw>& v5, const CompactCoordinatesXYZFF<vbw>& v6,
        const CompactCoordinatesXYZFF<vbw>& v7, const CompactCoordinatesXYZFF<vbw>& v8) const {
        return CompactCoordinatesXYZFF<vbw>::evaluate_rounded_sse(v1, v2, v3, v4, v5, v6, v7, v8);
    }
    #endif

    #if defined __AVX__
    xyzff::OctoEvaluatedResultRounded avx8(
        const CompactCoordinatesXYZFF<vbw>& v1, const CompactCoordinatesXYZFF<vbw>& v2,
        const CompactCoordinatesXYZFF<vbw>& v3, const CompactCoordinatesXYZFF<vbw>& v4,
        const CompactCoordinatesXYZFF<vbw>& v5, const CompactCoordinatesXYZFF<vbw>& v6,
        const CompactCoordinatesXYZFF<vbw>& v7, const CompactCoordinatesXYZFF<vbw>& v8) const {
        return CompactCoordinatesXYZFF<vbw>::evaluate_rounded_avx(v1, v2, v3, v4, v5, v6, v7, v8);
    }
    #endif
};

namespace {
    constexpr size_t N = 1024; // must be a multiple of 8

    template<bool vbw>
    std::vector<XYZWBench<vbw>> make_xyzw(unsigned seed = 42) {
        std::mt19937 gen(seed);
        std::uniform_real_distribution<float> pos(-100.f, 100.f);
        std::uniform_real_distribution<float> w(0.1f, 2.0f);
        std::vector<XYZWBench<vbw>> v;
        v.reserve(N);
        for (size_t i = 0; i < N; ++i) {
            v.emplace_back(Vector3<float>{pos(gen), pos(gen), pos(gen)}, w(gen));
        }
        return v;
    }

    template<bool vbw>
    std::vector<XYZFFBench<vbw>> make_xyzff(unsigned seed = 42) {
        std::mt19937 gen(seed);
        std::uniform_real_distribution<float> pos(-100.f, 100.f);
        std::uniform_int_distribution<int> ff(0, static_cast<int>(form_factor::form_factor_t::COUNT) - 1);
        std::vector<XYZFFBench<vbw>> v;
        v.reserve(N);
        for (size_t i = 0; i < N; ++i) {
            v.emplace_back(Vector3<float>{pos(gen), pos(gen), pos(gen)}, ff(gen));
        }
        return v;
    }
}

TEST_CASE("SIMD benchmark: CompactCoordinatesXYZW evaluate_rounded", "[bench]") {
    auto atoms = make_xyzw<false>();
    const auto& ref = atoms[0];

    BENCHMARK_ADVANCED("scalar (4 at a time)")(Catch::Benchmark::Chronometer meter) {
        int32_t sink = 0;
        meter.measure([&] {
            for (size_t j = 8; j + 4 <= N; j += 4) {
                auto r = ref.scalar4(atoms[j], atoms[j+1], atoms[j+2], atoms[j+3]);
                sink += r.distances[0];
            }
            return sink;
        });
    };

    BENCHMARK_ADVANCED("scalar (8 at a time)")(Catch::Benchmark::Chronometer meter) {
        int32_t sink = 0;
        meter.measure([&] {
            for (size_t j = 8; j + 8 <= N; j += 8) {
                auto r = ref.scalar8(atoms[j], atoms[j+1], atoms[j+2], atoms[j+3], atoms[j+4], atoms[j+5], atoms[j+6], atoms[j+7]);
                sink += r.distances[0];
            }
            return sink;
        });
    };

    #if defined __SSE2__
    BENCHMARK_ADVANCED("SSE (4 at a time)")(Catch::Benchmark::Chronometer meter) {
        int32_t sink = 0;
        meter.measure([&] {
            for (size_t j = 8; j + 4 <= N; j += 4) {
                auto r = ref.sse4(atoms[j], atoms[j+1], atoms[j+2], atoms[j+3]);
                sink += r.distances[0];
            }
            return sink;
        });
    };

    BENCHMARK_ADVANCED("SSE (8 at a time)")(Catch::Benchmark::Chronometer meter) {
        int32_t sink = 0;
        meter.measure([&] {
            for (size_t j = 8; j + 8 <= N; j += 8) {
                auto r = ref.sse8(atoms[j], atoms[j+1], atoms[j+2], atoms[j+3], atoms[j+4], atoms[j+5], atoms[j+6], atoms[j+7]);
                sink += r.distances[0];
            }
            return sink;
        });
    };
    #endif

    #if defined __AVX__
    BENCHMARK_ADVANCED("AVX (8 at a time)")(Catch::Benchmark::Chronometer meter) {
        int32_t sink = 0;
        meter.measure([&] {
            for (size_t j = 8; j + 8 <= N; j += 8) {
                auto r = ref.avx8(atoms[j], atoms[j+1], atoms[j+2], atoms[j+3], atoms[j+4], atoms[j+5], atoms[j+6], atoms[j+7]);
                sink += r.distances[0];
            }
            return sink;
        });
    };
    #endif
}

TEST_CASE("SIMD benchmark: CompactCoordinatesXYZFF evaluate_rounded", "[bench]") {
    auto atoms = make_xyzff<false>();
    const auto& ref = atoms[0];

    BENCHMARK_ADVANCED("scalar (4 at a time)")(Catch::Benchmark::Chronometer meter) {
        int32_t sink = 0;
        meter.measure([&] {
            for (size_t j = 8; j + 4 <= N; j += 4) {
                auto r = ref.scalar4(atoms[j], atoms[j+1], atoms[j+2], atoms[j+3]);
                sink += r.distances[0];
            }
            return sink;
        });
    };

    BENCHMARK_ADVANCED("scalar (8 at a time)")(Catch::Benchmark::Chronometer meter) {
        int32_t sink = 0;
        meter.measure([&] {
            for (size_t j = 8; j + 8 <= N; j += 8) {
                auto r = ref.scalar8(atoms[j], atoms[j+1], atoms[j+2], atoms[j+3], atoms[j+4], atoms[j+5], atoms[j+6], atoms[j+7]);
                sink += r.distances[0];
            }
            return sink;
        });
    };

    #if defined __SSE2__
    BENCHMARK_ADVANCED("SSE (4 at a time)")(Catch::Benchmark::Chronometer meter) {
        int32_t sink = 0;
        meter.measure([&] {
            for (size_t j = 8; j + 4 <= N; j += 4) {
                auto r = ref.sse4(atoms[j], atoms[j+1], atoms[j+2], atoms[j+3]);
                sink += r.distances[0];
            }
            return sink;
        });
    };

    BENCHMARK_ADVANCED("SSE (8 at a time)")(Catch::Benchmark::Chronometer meter) {
        int32_t sink = 0;
        meter.measure([&] {
            for (size_t j = 8; j + 8 <= N; j += 8) {
                auto r = ref.sse8(atoms[j], atoms[j+1], atoms[j+2], atoms[j+3], atoms[j+4], atoms[j+5], atoms[j+6], atoms[j+7]);
                sink += r.distances[0];
            }
            return sink;
        });
    };
    #endif

    #if defined __AVX__
    BENCHMARK_ADVANCED("AVX (8 at a time)")(Catch::Benchmark::Chronometer meter) {
        int32_t sink = 0;
        meter.measure([&] {
            for (size_t j = 8; j + 8 <= N; j += 8) {
                auto r = ref.avx8(atoms[j], atoms[j+1], atoms[j+2], atoms[j+3], atoms[j+4], atoms[j+5], atoms[j+6], atoms[j+7]);
                sink += r.distances[0];
            }
            return sink;
        });
    };
    #endif
}
