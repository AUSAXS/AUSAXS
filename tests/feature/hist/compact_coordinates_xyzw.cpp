#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <hist/detail/data/CompactCoordinatesXYZW.h>
#include <constants/Constants.h>
#include <math/Vector3.h>

#include <algorithm>
#include <array>
#include <functional>
#include <numeric>
#include <span>
#include <vector>

using namespace ausaxs;
using namespace hist::detail;
using namespace hist::detail::xyzw;

TEST_CASE("CompactCoordinatesXYZW<vbw>::CompactCoordinatesData") {
    SECTION("Vector3<double>, float") {
        CompactCoordinatesXYZW<false> data(Vector3<double>(1, 2, 3), 4);
        CHECK(data.value.pos.x() == 1);
        CHECK(data.value.pos.y() == 2);
        CHECK(data.value.pos.z() == 3);
        CHECK(data.value.w == 4);
    }
}

template<bool vbw>
struct DebugData : CompactCoordinatesXYZW<vbw> {
    using CC = CompactCoordinatesXYZW<vbw>;
    using CC::CC;

    EvaluatedResult evaluate_scalar(const CC& other) const {return CC::evaluate_scalar(other);}
    QuadEvaluatedResult evaluate_4_scalar(std::span<const CC, 4> others) const {return CC::evaluate_4_scalar(others);}
    OctoEvaluatedResult evaluate_8_scalar(std::span<const CC, 8> others) const {return CC::evaluate_8(others);}

    EvaluatedResultRounded evaluate_rounded_scalar(const CC& other) const {return CC::evaluate_rounded_scalar(other);}
    QuadEvaluatedResultRounded evaluate_rounded_4_scalar(std::span<const CC, 4> others) const {return CC::evaluate_rounded_4_scalar(others);}
    OctoEvaluatedResultRounded evaluate_rounded_8_scalar(std::span<const CC, 8> others) const {return CC::evaluate_rounded_8(others);}

    #if defined AUSAXS_USE_SSE2
        QuadEvaluatedResult evaluate_4_sse(std::span<const CC, 4> others) const {return CC::evaluate_4_sse(others);}
        QuadEvaluatedResultRounded evaluate_rounded_4_sse(std::span<const CC, 4> others) const {return CC::evaluate_rounded_4_sse(others);}
    #endif

    #if defined AUSAXS_USE_AVX2
        OctoEvaluatedResult evaluate_8_avx(std::span<const CC, 8> others) const {return CC::evaluate_8_avx(others);}
        OctoEvaluatedResultRounded evaluate_rounded_8_avx(std::span<const CC, 8> others) const {return CC::evaluate_rounded_8_avx(others);}
    #endif

    #if defined AUSAXS_USE_AVX512
        HexaEvaluatedResult evaluate_16_avx512(std::span<const CC, 16> others) const {return CC::evaluate_16_avx512(others);}
        HexaEvaluatedResultRounded evaluate_rounded_16_avx512(std::span<const CC, 16> others) const {return CC::evaluate_rounded_16_avx512(others);}
    #endif
};

template<bool vbw>
using CC = CompactCoordinatesXYZW<vbw>;

// SIMD backends may reorder output elements; sort by distance and compare as sets
template<std::size_t N>
void check_unordered(
    const std::array<float, N>& actual_dist,
    const std::array<float, N>& actual_wt,
    std::vector<std::pair<double, float>> expected,
    double tol)
{
    std::sort(expected.begin(), expected.end());
    std::array<std::size_t, N> idx;
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](auto a, auto b) { return actual_dist[a] < actual_dist[b]; });
    for (std::size_t k = 0; k < N; ++k) {
        CHECK_THAT(static_cast<double>(actual_dist[idx[k]]), Catch::Matchers::WithinAbs(expected[k].first, tol));
        CHECK(actual_wt[idx[k]] == expected[k].second);
    }
}

template<std::size_t N>
void check_unordered_rounded(
    const std::array<int32_t, N>& actual_dist,
    const std::array<float, N>& actual_wt,
    std::vector<std::pair<int32_t, float>> expected)
{
    std::sort(expected.begin(), expected.end());
    std::array<std::size_t, N> idx;
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](auto a, auto b) { return actual_dist[a] < actual_dist[b]; });
    for (std::size_t k = 0; k < N; ++k) {
        CHECK(actual_dist[idx[k]] == expected[k].first);
        CHECK(actual_wt[idx[k]] == expected[k].second);
    }
}

template<bool vbw>
void single_tests(std::function<EvaluatedResult(const DebugData<vbw>&, const DebugData<vbw>&)> evaluate) {
    SECTION("single distance") {
        DebugData<vbw> data1(Vector3<double>{1, 1, 1}, 2);
        DebugData<vbw> data2(Vector3<double>{2, 1, 1}, 4);
        auto result = evaluate(data1, data2);
        CHECK(result.distance == 1);
        CHECK(result.weight == 8);

        DebugData<vbw> data3(Vector3<double>{2, 2, 2}, 8);
        result = evaluate(data1, data3);
        CHECK_THAT(result.distance, Catch::Matchers::WithinAbs(std::sqrt(3), 1e-6));
        CHECK(result.weight == 16);
    }
}

template<bool vbw>
void single_tests_rounded(std::function<EvaluatedResultRounded(const DebugData<vbw>&, const DebugData<vbw>&)> evaluate) {
    SECTION("single distance") {
        double width = constants::axes::d_axis.width();
        DebugData<vbw> data1(Vector3<double>{1, 1, 1}, 2);
        DebugData<vbw> data2(Vector3<double>{2, 1, 1}, 4);
        auto result = evaluate(data1, data2);
        CHECK(result.distance == std::round(1./width));
        CHECK(result.weight == 8);

        DebugData<vbw> data3(Vector3<double>{2, 2, 2}, 8);
        result = evaluate(data1, data3);
        CHECK(result.distance == std::round(std::sqrt(3)/width));
        CHECK(result.weight == 16);
    }
}

template<bool vbw>
void quad_tests(std::function<QuadEvaluatedResult(const DebugData<vbw>&, const std::array<CC<vbw>, 4>&)> evaluate) {
    SECTION("four distances") {
        DebugData<vbw> data(Vector3<double>{1, 1, 1}, 2);
        std::array<CC<vbw>, 4> others = {
            CC<vbw>(Vector3<double>{2, 1, 1}, 4),
            CC<vbw>(Vector3<double>{2, 2, 2}, 8),
            CC<vbw>(Vector3<double>{3, 3, 3}, 16),
            CC<vbw>(Vector3<double>{4, 4, 4}, 3)
        };
        auto result = evaluate(data, others);
        check_unordered<4>(result.distances, result.weights, {
            {1.0, 8.0f}, {std::sqrt(3.0), 16.0f}, {std::sqrt(12.0), 32.0f}, {std::sqrt(27.0), 6.0f}
        }, 1e-6);
    }
}

template<bool vbw>
void quad_tests_rounded(std::function<QuadEvaluatedResultRounded(const DebugData<vbw>&, const std::array<CC<vbw>, 4>&)> evaluate) {
    SECTION("four distances") {
        double width = constants::axes::d_axis.width();
        DebugData<vbw> data(Vector3<double>{1, 1, 1}, 2);
        std::array<CC<vbw>, 4> others = {
            CC<vbw>(Vector3<double>{2, 1, 1}, 4),
            CC<vbw>(Vector3<double>{2, 2, 2}, 8),
            CC<vbw>(Vector3<double>{3, 3, 3}, 16),
            CC<vbw>(Vector3<double>{4, 4, 4}, 3)
        };
        auto result = evaluate(data, others);
        check_unordered_rounded<4>(result.distances, result.weights, {
            {static_cast<int32_t>(std::round(1.0/width)), 8.0f},
            {static_cast<int32_t>(std::round(std::sqrt(3.0)/width)), 16.0f},
            {static_cast<int32_t>(std::round(std::sqrt(12.0)/width)), 32.0f},
            {static_cast<int32_t>(std::round(std::sqrt(27.0)/width)), 6.0f}
        });
    }
}

template<bool vbw>
void octo_tests(std::function<OctoEvaluatedResult(const DebugData<vbw>&, const std::array<CC<vbw>, 8>&)> evaluate) {
    SECTION("eight distances") {
        DebugData<vbw> data(Vector3<double>{1, 1, 1}, 2);
        std::array<CC<vbw>, 8> others = {
            CC<vbw>(Vector3<double>{2, 1, 1}, 4),
            CC<vbw>(Vector3<double>{2, 2, 2}, 8),
            CC<vbw>(Vector3<double>{3, 3, 3}, 16),
            CC<vbw>(Vector3<double>{4, 4, 4}, 32),
            CC<vbw>(Vector3<double>{5, 5, 5}, 64),
            CC<vbw>(Vector3<double>{6, 6, 6}, 128),
            CC<vbw>(Vector3<double>{7, 7, 7}, 15),
            CC<vbw>(Vector3<double>{8, 8, 8}, 5)
        };
        auto result = evaluate(data, others);
        check_unordered<8>(result.distances, result.weights, {
            {1.0, 8.0f}, {std::sqrt(3.0), 16.0f}, {std::sqrt(12.0), 32.0f}, {std::sqrt(27.0), 64.0f},
            {std::sqrt(48.0), 128.0f}, {std::sqrt(75.0), 256.0f}, {std::sqrt(108.0), 30.0f}, {std::sqrt(147.0), 10.0f}
        }, 1e-5);
    }
}

template<bool vbw>
void octo_tests_rounded(std::function<OctoEvaluatedResultRounded(const DebugData<vbw>&, const std::array<CC<vbw>, 8>&)> evaluate) {
    SECTION("eight distances") {
        double width = constants::axes::d_axis.width();
        DebugData<vbw> data(Vector3<double>{1, 1, 1}, 2);
        std::array<CC<vbw>, 8> others = {
            CC<vbw>(Vector3<double>{2, 1, 1}, 4),
            CC<vbw>(Vector3<double>{2, 2, 2}, 8),
            CC<vbw>(Vector3<double>{3, 3, 3}, 16),
            CC<vbw>(Vector3<double>{4, 4, 4}, 32),
            CC<vbw>(Vector3<double>{5, 5, 5}, 64),
            CC<vbw>(Vector3<double>{6, 6, 6}, 128),
            CC<vbw>(Vector3<double>{7, 7, 7}, 15),
            CC<vbw>(Vector3<double>{8, 8, 8}, 5)
        };
        auto result = evaluate(data, others);
        check_unordered_rounded<8>(result.distances, result.weights, {
            {static_cast<int32_t>(std::round(1.0/width)), 8.0f},
            {static_cast<int32_t>(std::round(std::sqrt(3.0)/width)), 16.0f},
            {static_cast<int32_t>(std::round(std::sqrt(12.0)/width)), 32.0f},
            {static_cast<int32_t>(std::round(std::sqrt(27.0)/width)), 64.0f},
            {static_cast<int32_t>(std::round(std::sqrt(48.0)/width)), 128.0f},
            {static_cast<int32_t>(std::round(std::sqrt(75.0)/width)), 256.0f},
            {static_cast<int32_t>(std::round(std::sqrt(108.0)/width)), 30.0f},
            {static_cast<int32_t>(std::round(std::sqrt(147.0)/width)), 10.0f}
        });
    }
}

template<bool vbw>
void hexa_tests(std::function<HexaEvaluatedResult(const DebugData<vbw>&, const std::array<CC<vbw>, 16>&)> evaluate) {
    SECTION("sixteen distances") {
        DebugData<vbw> data(Vector3<double>{1, 1, 1}, 2);
        std::array<CC<vbw>, 16> others = {
            CC<vbw>(Vector3<double>{2, 1, 1}, 4),
            CC<vbw>(Vector3<double>{2, 2, 2}, 8),
            CC<vbw>(Vector3<double>{3, 3, 3}, 16),
            CC<vbw>(Vector3<double>{4, 4, 4}, 32),
            CC<vbw>(Vector3<double>{5, 5, 5}, 64),
            CC<vbw>(Vector3<double>{6, 6, 6}, 128),
            CC<vbw>(Vector3<double>{7, 7, 7}, 15),
            CC<vbw>(Vector3<double>{8, 8, 8}, 5),
            CC<vbw>(Vector3<double>{9, 9, 9}, 7),
            CC<vbw>(Vector3<double>{10, 10, 10}, 11),
            CC<vbw>(Vector3<double>{11, 11, 11}, 13),
            CC<vbw>(Vector3<double>{12, 12, 12}, 17),
            CC<vbw>(Vector3<double>{13, 13, 13}, 19),
            CC<vbw>(Vector3<double>{14, 14, 14}, 23),
            CC<vbw>(Vector3<double>{15, 15, 15}, 29),
            CC<vbw>(Vector3<double>{16, 16, 16}, 31)
        };
        auto result = evaluate(data, others);
        check_unordered<16>(result.distances, result.weights, {
            {1.0, 8.0f}, {std::sqrt(3.0), 16.0f}, {std::sqrt(12.0), 32.0f}, {std::sqrt(27.0), 64.0f},
            {std::sqrt(48.0), 128.0f}, {std::sqrt(75.0), 256.0f}, {std::sqrt(108.0), 30.0f}, {std::sqrt(147.0), 10.0f},
            {std::sqrt(192.0), 14.0f}, {std::sqrt(243.0), 22.0f}, {std::sqrt(300.0), 26.0f}, {std::sqrt(363.0), 34.0f},
            {std::sqrt(432.0), 38.0f}, {std::sqrt(507.0), 46.0f}, {std::sqrt(588.0), 58.0f}, {std::sqrt(675.0), 62.0f}
        }, 1e-3);
    }
}

template<bool vbw>
void hexa_tests_rounded(std::function<HexaEvaluatedResultRounded(const DebugData<vbw>&, const std::array<CC<vbw>, 16>&)> evaluate) {
    SECTION("sixteen distances") {
        double width = constants::axes::d_axis.width();
        DebugData<vbw> data(Vector3<double>{1, 1, 1}, 2);
        std::array<CC<vbw>, 16> others = {
            CC<vbw>(Vector3<double>{2, 1, 1}, 4),
            CC<vbw>(Vector3<double>{2, 2, 2}, 8),
            CC<vbw>(Vector3<double>{3, 3, 3}, 16),
            CC<vbw>(Vector3<double>{4, 4, 4}, 32),
            CC<vbw>(Vector3<double>{5, 5, 5}, 64),
            CC<vbw>(Vector3<double>{6, 6, 6}, 128),
            CC<vbw>(Vector3<double>{7, 7, 7}, 15),
            CC<vbw>(Vector3<double>{8, 8, 8}, 5),
            CC<vbw>(Vector3<double>{9, 9, 9}, 7),
            CC<vbw>(Vector3<double>{10, 10, 10}, 11),
            CC<vbw>(Vector3<double>{11, 11, 11}, 13),
            CC<vbw>(Vector3<double>{12, 12, 12}, 17),
            CC<vbw>(Vector3<double>{13, 13, 13}, 19),
            CC<vbw>(Vector3<double>{14, 14, 14}, 23),
            CC<vbw>(Vector3<double>{15, 15, 15}, 29),
            CC<vbw>(Vector3<double>{16, 16, 16}, 31)
        };
        auto result = evaluate(data, others);
        check_unordered_rounded<16>(result.distances, result.weights, {
            {static_cast<int32_t>(std::round(1.0/width)), 8.0f},
            {static_cast<int32_t>(std::round(std::sqrt(3.0)/width)), 16.0f},
            {static_cast<int32_t>(std::round(std::sqrt(12.0)/width)), 32.0f},
            {static_cast<int32_t>(std::round(std::sqrt(27.0)/width)), 64.0f},
            {static_cast<int32_t>(std::round(std::sqrt(48.0)/width)), 128.0f},
            {static_cast<int32_t>(std::round(std::sqrt(75.0)/width)), 256.0f},
            {static_cast<int32_t>(std::round(std::sqrt(108.0)/width)), 30.0f},
            {static_cast<int32_t>(std::round(std::sqrt(147.0)/width)), 10.0f},
            {static_cast<int32_t>(std::round(std::sqrt(192.0)/width)), 14.0f},
            {static_cast<int32_t>(std::round(std::sqrt(243.0)/width)), 22.0f},
            {static_cast<int32_t>(std::round(std::sqrt(300.0)/width)), 26.0f},
            {static_cast<int32_t>(std::round(std::sqrt(363.0)/width)), 34.0f},
            {static_cast<int32_t>(std::round(std::sqrt(432.0)/width)), 38.0f},
            {static_cast<int32_t>(std::round(std::sqrt(507.0)/width)), 46.0f},
            {static_cast<int32_t>(std::round(std::sqrt(588.0)/width)), 58.0f},
            {static_cast<int32_t>(std::round(std::sqrt(675.0)/width)), 62.0f}
        });
    }
}

template<bool vbw>
void run_tests() {
    SECTION("scalar") {
        single_tests<vbw>([](const DebugData<vbw>& d1, const DebugData<vbw>& d2) { return d1.evaluate_scalar(d2); });
        quad_tests<vbw>([](const DebugData<vbw>& d, const std::array<CC<vbw>, 4>& o) { return d.evaluate_4_scalar(std::span<const CC<vbw>, 4>(o)); });
        octo_tests<vbw>([](const DebugData<vbw>& d, const std::array<CC<vbw>, 8>& o) { return d.evaluate_8_scalar(std::span<const CC<vbw>, 8>(o)); });

        single_tests_rounded<vbw>([](const DebugData<vbw>& d1, const DebugData<vbw>& d2) { return d1.evaluate_rounded_scalar(d2); });
        quad_tests_rounded<vbw>([](const DebugData<vbw>& d, const std::array<CC<vbw>, 4>& o) { return d.evaluate_rounded_4_scalar(std::span<const CC<vbw>, 4>(o)); });
        octo_tests_rounded<vbw>([](const DebugData<vbw>& d, const std::array<CC<vbw>, 8>& o) { return d.evaluate_rounded_8_scalar(std::span<const CC<vbw>, 8>(o)); });
    }

    #if defined AUSAXS_USE_SSE2
        SECTION("sse") {
            quad_tests<vbw>([](const DebugData<vbw>& d, const std::array<CC<vbw>, 4>& o) { return d.evaluate_4_sse(std::span<const CC<vbw>, 4>(o)); });
            quad_tests_rounded<vbw>([](const DebugData<vbw>& d, const std::array<CC<vbw>, 4>& o) { return d.evaluate_rounded_4_sse(std::span<const CC<vbw>, 4>(o)); });
        }
    #endif

    #if defined AUSAXS_USE_AVX2
        SECTION("avx") {
            octo_tests<vbw>([](const DebugData<vbw>& d, const std::array<CC<vbw>, 8>& o) { return d.evaluate_8_avx(std::span<const CC<vbw>, 8>(o)); });
            octo_tests_rounded<vbw>([](const DebugData<vbw>& d, const std::array<CC<vbw>, 8>& o) { return d.evaluate_rounded_8_avx(std::span<const CC<vbw>, 8>(o)); });
        }
    #endif

    #if defined AUSAXS_USE_AVX512
        SECTION("avx512") {
            hexa_tests<vbw>([](const DebugData<vbw>& d, const std::array<CC<vbw>, 16>& o) { return d.evaluate_16_avx512(std::span<const CC<vbw>, 16>(o)); });
            hexa_tests_rounded<vbw>([](const DebugData<vbw>& d, const std::array<CC<vbw>, 16>& o) { return d.evaluate_rounded_16_avx512(std::span<const CC<vbw>, 16>(o)); });
        }
    #endif

    SECTION("dispatch") {
        hexa_tests<vbw>([](const DebugData<vbw>& d, const std::array<CC<vbw>, 16>& o) { return d.evaluate_16(std::span<const CC<vbw>, 16>(o)); });
        hexa_tests_rounded<vbw>([](const DebugData<vbw>& d, const std::array<CC<vbw>, 16>& o) { return d.evaluate_rounded_16(std::span<const CC<vbw>, 16>(o)); });
    }
}

TEST_CASE("CompactCoordinatesXYZW<vbw>::evaluate") {
    SECTION("variable bin width") {
        run_tests<true>();
    }
    SECTION("fixed bin width") {
        run_tests<false>();
    }
}