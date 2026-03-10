#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <hist/detail/data/CompactCoordinatesXYZFF.h>
#include <constants/Constants.h>
#include <math/Vector3.h>
#include <form_factor/FormFactorType.h>

#include <algorithm>
#include <array>
#include <functional>
#include <numeric>
#include <span>
#include <vector>

using namespace ausaxs;
using namespace hist::detail;
using namespace hist::detail::xyzff;

TEST_CASE("CompactCoordinatesXYZFF<vbw>::CompactCoordinatesData") {
    SECTION("Vector3<double>, int32_t") {
        CompactCoordinatesXYZFF<false> data(Vector3<double>(1, 2, 3), 4);
        CHECK(data.value.pos.x() == 1);
        CHECK(data.value.pos.y() == 2);
        CHECK(data.value.pos.z() == 3);
        CHECK(data.value.ff == 4);
    }
}

template<bool vbw>
struct DebugData : CompactCoordinatesXYZFF<vbw> {
    using CC = CompactCoordinatesXYZFF<vbw>;
    using CC::CC;

    EvaluatedResult evaluate_scalar(const CC& other) const {return CC::evaluate_scalar(other);}
    QuadEvaluatedResult evaluate_4_scalar(std::span<const CC, 4> others) const {return CC::evaluate_4_scalar(others);}
    OctoEvaluatedResult evaluate_8_scalar(std::span<const CC, 8> others) const {return CC::evaluate_8(others);}

    EvaluatedResultRounded evaluate_rounded_scalar(const CC& other) const {return CC::evaluate_rounded_scalar(other);}
    QuadEvaluatedResultRounded evaluate_rounded_4_scalar(std::span<const CC, 4> others) const {return CC::evaluate_rounded_4_scalar(others);}

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
using CC = CompactCoordinatesXYZFF<vbw>;

// SIMD backends may reorder output elements; sort by distance and compare as sets
template<std::size_t N>
void check_unordered(
    const std::array<float, N>& actual_dist,
    const std::array<int32_t, N>& actual_ff,
    std::vector<std::pair<double, int32_t>> expected,
    double tol)
{
    std::sort(expected.begin(), expected.end());
    std::array<std::size_t, N> idx;
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](auto a, auto b) { return actual_dist[a] < actual_dist[b]; });
    for (std::size_t k = 0; k < N; ++k) {
        CHECK_THAT(static_cast<double>(actual_dist[idx[k]]), Catch::Matchers::WithinAbs(expected[k].first, tol));
        CHECK(actual_ff[idx[k]] == expected[k].second);
    }
}

template<std::size_t N>
void check_unordered_rounded(
    const std::array<int32_t, N>& actual_dist,
    const std::array<int32_t, N>& actual_ff,
    std::vector<std::pair<int32_t, int32_t>> expected)
{
    std::sort(expected.begin(), expected.end());
    std::array<std::size_t, N> idx;
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](auto a, auto b) { return actual_dist[a] < actual_dist[b]; });
    for (std::size_t k = 0; k < N; ++k) {
        CHECK(actual_dist[idx[k]] == expected[k].first);
        CHECK(actual_ff[idx[k]] == expected[k].second);
    }
}

template<bool vbw>
void single_tests(std::function<EvaluatedResult(const DebugData<vbw>&, const DebugData<vbw>&)> evaluate) {
    SECTION("single distance") {
        DebugData<vbw> data1(Vector3<double>{1, 1, 1}, 2);
        DebugData<vbw> data2(Vector3<double>{2, 1, 1}, 4);
        auto result = evaluate(data1, data2);
        CHECK(result.distance == 1);
        CHECK(result.ff_bin == ff_bin_index<false>(2, 4));

        DebugData<vbw> data3(Vector3<double>{2, 2, 2}, 8);
        result = evaluate(data1, data3);
        CHECK_THAT(result.distance, Catch::Matchers::WithinAbs(std::sqrt(3), 1e-6));
        CHECK(result.ff_bin == ff_bin_index<false>(2, 8));
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
        CHECK(result.ff_bin == ff_bin_index<false>(2, 4));

        DebugData<vbw> data3(Vector3<double>{2, 2, 2}, 8);
        result = evaluate(data1, data3);
        CHECK(result.distance == std::round(std::sqrt(3)/width));
        CHECK(result.ff_bin == ff_bin_index<false>(2, 8));
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
        check_unordered<4>(result.distances, result.ff_bins, {
            {1.0, ff_bin_index<false>(2, 4)}, {std::sqrt(3.0), ff_bin_index<false>(2, 8)},
            {std::sqrt(12.0), ff_bin_index<false>(2, 16)}, {std::sqrt(27.0), ff_bin_index<false>(2, 3)}
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
        check_unordered_rounded<4>(result.distances, result.ff_bins, {
            {static_cast<int32_t>(std::round(1.0/width)), ff_bin_index<false>(2, 4)},
            {static_cast<int32_t>(std::round(std::sqrt(3.0)/width)), ff_bin_index<false>(2, 8)},
            {static_cast<int32_t>(std::round(std::sqrt(12.0)/width)), ff_bin_index<false>(2, 16)},
            {static_cast<int32_t>(std::round(std::sqrt(27.0)/width)), ff_bin_index<false>(2, 3)}
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
        check_unordered<8>(result.distances, result.ff_bins, {
            {1.0, ff_bin_index<false>(2, 4)}, {std::sqrt(3.0), ff_bin_index<false>(2, 8)},
            {std::sqrt(12.0), ff_bin_index<false>(2, 16)}, {std::sqrt(27.0), ff_bin_index<false>(2, 32)},
            {std::sqrt(48.0), ff_bin_index<false>(2, 64)}, {std::sqrt(75.0), ff_bin_index<false>(2, 128)},
            {std::sqrt(108.0), ff_bin_index<false>(2, 15)}, {std::sqrt(147.0), ff_bin_index<false>(2, 5)}
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
        check_unordered_rounded<8>(result.distances, result.ff_bins, {
            {static_cast<int32_t>(std::round(1.0/width)), ff_bin_index<false>(2, 4)},
            {static_cast<int32_t>(std::round(std::sqrt(3.0)/width)), ff_bin_index<false>(2, 8)},
            {static_cast<int32_t>(std::round(std::sqrt(12.0)/width)), ff_bin_index<false>(2, 16)},
            {static_cast<int32_t>(std::round(std::sqrt(27.0)/width)), ff_bin_index<false>(2, 32)},
            {static_cast<int32_t>(std::round(std::sqrt(48.0)/width)), ff_bin_index<false>(2, 64)},
            {static_cast<int32_t>(std::round(std::sqrt(75.0)/width)), ff_bin_index<false>(2, 128)},
            {static_cast<int32_t>(std::round(std::sqrt(108.0)/width)), ff_bin_index<false>(2, 15)},
            {static_cast<int32_t>(std::round(std::sqrt(147.0)/width)), ff_bin_index<false>(2, 5)}
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
            CC<vbw>(Vector3<double>{3, 3, 3}, 3),
            CC<vbw>(Vector3<double>{4, 4, 4}, 5),
            CC<vbw>(Vector3<double>{5, 5, 5}, 6),
            CC<vbw>(Vector3<double>{6, 6, 6}, 7),
            CC<vbw>(Vector3<double>{7, 7, 7}, 1),
            CC<vbw>(Vector3<double>{8, 8, 8}, 9),
            CC<vbw>(Vector3<double>{9, 9, 9}, 10),
            CC<vbw>(Vector3<double>{10, 10, 10}, 11),
            CC<vbw>(Vector3<double>{11, 11, 11}, 12),
            CC<vbw>(Vector3<double>{12, 12, 12}, 13),
            CC<vbw>(Vector3<double>{13, 13, 13}, 14),
            CC<vbw>(Vector3<double>{14, 14, 14}, 15),
            CC<vbw>(Vector3<double>{15, 15, 15}, 16),
            CC<vbw>(Vector3<double>{16, 16, 16}, 0)
        };
        auto result = evaluate(data, others);
        check_unordered<16>(result.distances, result.ff_bins, {
            {1.0, ff_bin_index<false>(2, 4)}, {std::sqrt(3.0), ff_bin_index<false>(2, 8)},
            {std::sqrt(12.0), ff_bin_index<false>(2, 3)}, {std::sqrt(27.0), ff_bin_index<false>(2, 5)},
            {std::sqrt(48.0), ff_bin_index<false>(2, 6)}, {std::sqrt(75.0), ff_bin_index<false>(2, 7)},
            {std::sqrt(108.0), ff_bin_index<false>(2, 1)}, {std::sqrt(147.0), ff_bin_index<false>(2, 9)},
            {std::sqrt(192.0), ff_bin_index<false>(2, 10)}, {std::sqrt(243.0), ff_bin_index<false>(2, 11)},
            {std::sqrt(300.0), ff_bin_index<false>(2, 12)}, {std::sqrt(363.0), ff_bin_index<false>(2, 13)},
            {std::sqrt(432.0), ff_bin_index<false>(2, 14)}, {std::sqrt(507.0), ff_bin_index<false>(2, 15)},
            {std::sqrt(588.0), ff_bin_index<false>(2, 16)}, {std::sqrt(675.0), ff_bin_index<false>(2, 0)}
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
            CC<vbw>(Vector3<double>{3, 3, 3}, 3),
            CC<vbw>(Vector3<double>{4, 4, 4}, 5),
            CC<vbw>(Vector3<double>{5, 5, 5}, 6),
            CC<vbw>(Vector3<double>{6, 6, 6}, 7),
            CC<vbw>(Vector3<double>{7, 7, 7}, 1),
            CC<vbw>(Vector3<double>{8, 8, 8}, 9),
            CC<vbw>(Vector3<double>{9, 9, 9}, 10),
            CC<vbw>(Vector3<double>{10, 10, 10}, 11),
            CC<vbw>(Vector3<double>{11, 11, 11}, 12),
            CC<vbw>(Vector3<double>{12, 12, 12}, 13),
            CC<vbw>(Vector3<double>{13, 13, 13}, 14),
            CC<vbw>(Vector3<double>{14, 14, 14}, 15),
            CC<vbw>(Vector3<double>{15, 15, 15}, 16),
            CC<vbw>(Vector3<double>{16, 16, 16}, 0)
        };
        auto result = evaluate(data, others);
        check_unordered_rounded<16>(result.distances, result.ff_bins, {
            {static_cast<int32_t>(std::round(1.0/width)), ff_bin_index<false>(2, 4)},
            {static_cast<int32_t>(std::round(std::sqrt(3.0)/width)), ff_bin_index<false>(2, 8)},
            {static_cast<int32_t>(std::round(std::sqrt(12.0)/width)), ff_bin_index<false>(2, 3)},
            {static_cast<int32_t>(std::round(std::sqrt(27.0)/width)), ff_bin_index<false>(2, 5)},
            {static_cast<int32_t>(std::round(std::sqrt(48.0)/width)), ff_bin_index<false>(2, 6)},
            {static_cast<int32_t>(std::round(std::sqrt(75.0)/width)), ff_bin_index<false>(2, 7)},
            {static_cast<int32_t>(std::round(std::sqrt(108.0)/width)), ff_bin_index<false>(2, 1)},
            {static_cast<int32_t>(std::round(std::sqrt(147.0)/width)), ff_bin_index<false>(2, 9)},
            {static_cast<int32_t>(std::round(std::sqrt(192.0)/width)), ff_bin_index<false>(2, 10)},
            {static_cast<int32_t>(std::round(std::sqrt(243.0)/width)), ff_bin_index<false>(2, 11)},
            {static_cast<int32_t>(std::round(std::sqrt(300.0)/width)), ff_bin_index<false>(2, 12)},
            {static_cast<int32_t>(std::round(std::sqrt(363.0)/width)), ff_bin_index<false>(2, 13)},
            {static_cast<int32_t>(std::round(std::sqrt(432.0)/width)), ff_bin_index<false>(2, 14)},
            {static_cast<int32_t>(std::round(std::sqrt(507.0)/width)), ff_bin_index<false>(2, 15)},
            {static_cast<int32_t>(std::round(std::sqrt(588.0)/width)), ff_bin_index<false>(2, 16)},
            {static_cast<int32_t>(std::round(std::sqrt(675.0)/width)), ff_bin_index<false>(2, 0)}
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

TEST_CASE("CompactCoordinatesXYZFF<vbw>::evaluate") {
    SECTION("variable bin width") {
        run_tests<true>();
    }
    SECTION("fixed bin width") {
        run_tests<false>();
    }
}
