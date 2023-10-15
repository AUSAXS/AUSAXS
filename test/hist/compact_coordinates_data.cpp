#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <hist/detail/CompactCoordinatesData.h>
#include <math/Vector3.h>

using namespace hist::detail;

TEST_CASE("CompactCoordinatesData::CompactCoordinatesData") {
    SECTION("default") {
        CompactCoordinatesData data;
        CHECK(data.value.x == 0);
        CHECK(data.value.y == 0);
        CHECK(data.value.z == 0);
        CHECK(data.value.w == 0);
    }

    SECTION("Vector3<double>, float") {
        CompactCoordinatesData data(Vector3<double>(1, 2, 3), 4);
        CHECK(data.value.x == 1);
        CHECK(data.value.y == 2);
        CHECK(data.value.z == 3);
        CHECK(data.value.w == 4);
    }
}

struct DebugData : CompactCoordinatesData {
    using CompactCoordinatesData::CompactCoordinatesData;
    EvaluatedResult evaluate_scalar(const CompactCoordinatesData& other) const {return CompactCoordinatesData::evaluate_scalar(other);}
    QuadEvaluatedResult evaluate_scalar(const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4) const {return CompactCoordinatesData::evaluate_scalar(v1, v2, v3, v4);}
    OctoEvaluatedResult evaluate_scalar(const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4, const CompactCoordinatesData& v5, const CompactCoordinatesData& v6, const CompactCoordinatesData& v7, const CompactCoordinatesData& v8) const {return CompactCoordinatesData::evaluate_scalar(v1, v2, v3, v4, v5, v6, v7, v8);}

    #if defined __SSE2__
        EvaluatedResult evaluate_sse(const CompactCoordinatesData& other) const {return CompactCoordinatesData::evaluate_sse(other);}
        QuadEvaluatedResult evaluate_sse(const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4) const {return CompactCoordinatesData::evaluate_sse(v1, v2, v3, v4);}
        OctoEvaluatedResult evaluate_sse(const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4, const CompactCoordinatesData& v5, const CompactCoordinatesData& v6, const CompactCoordinatesData& v7, const CompactCoordinatesData& v8) const {return CompactCoordinatesData::evaluate_sse(v1, v2, v3, v4, v5, v6, v7, v8);}
    #endif

    #if defined __AVX__
        EvaluatedResult evaluate_avx(const CompactCoordinatesData& other) const {return CompactCoordinatesData::evaluate_avx(other);}
        QuadEvaluatedResult evaluate_avx(const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4) const {return CompactCoordinatesData::evaluate_avx(v1, v2, v3, v4);}
        OctoEvaluatedResult evaluate_avx(const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4, const CompactCoordinatesData& v5, const CompactCoordinatesData& v6, const CompactCoordinatesData& v7, const CompactCoordinatesData& v8) const {return CompactCoordinatesData::evaluate_avx(v1, v2, v3, v4, v5, v6, v7, v8);}
    #endif
};

#include <functional>
void single_tests(std::function<EvaluatedResult(const DebugData&, const DebugData&)> evaluate) {
    DYNAMIC_SECTION("single distance, w = " << CompactCoordinatesData::inv_width) {
        double width = CompactCoordinatesData::inv_width;
        DebugData data1({1, 1, 1}, 2);
        DebugData data2({2, 1, 1}, 4);
        auto result = evaluate(data1, data2);
        CHECK(result.distance == std::round(1*width));
        CHECK(result.weight == 8);

        DebugData data3({2, 2, 2}, 8);
        result = evaluate(data1, data3);
        CHECK(result.distance == std::round(std::sqrt(3)*width));
        CHECK(result.weight == 16);
    }
}

void quad_tests(std::function<QuadEvaluatedResult(const DebugData&, const DebugData&, const DebugData&, const DebugData&, const DebugData&)> evaluate) {
    DYNAMIC_SECTION("four distances, w = " << CompactCoordinatesData::inv_width) {
        double width = CompactCoordinatesData::inv_width;
        DebugData data({1, 1, 1}, 2);
        DebugData data1({2, 1, 1}, 4);
        DebugData data2({2, 2, 2}, 8);
        DebugData data3({3, 3, 3}, 16);
        DebugData data4({4, 4, 4}, 3);
        auto result = evaluate(data, data1, data2, data3, data4);
        CHECK(result.distances.first  == std::round(1*width));
        CHECK(result.distances.second == std::round(std::sqrt(3)*width));
        CHECK(result.distances.third  == std::round(std::sqrt(12)*width));
        CHECK(result.distances.fourth == std::round(std::sqrt(27)*width));
        CHECK(result.weights.first == 8);
        CHECK(result.weights.second == 16);
        CHECK(result.weights.third == 32);
        CHECK(result.weights.fourth == 6);
    }
}

void octo_tests(std::function<OctoEvaluatedResult(const DebugData&, const DebugData&, const DebugData&, const DebugData&, const DebugData&, const DebugData&, const DebugData&, const DebugData&, const DebugData&)> evaluate) {
    DYNAMIC_SECTION("eight distances, w = " << CompactCoordinatesData::inv_width) {
        double width = CompactCoordinatesData::inv_width;
        DebugData data({1, 1, 1}, 2);
        DebugData data1({2, 1, 1}, 4);
        DebugData data2({2, 2, 2}, 8);
        DebugData data3({3, 3, 3}, 16);
        DebugData data4({4, 4, 4}, 32);
        DebugData data5({5, 5, 5}, 64);
        DebugData data6({6, 6, 6}, 128);
        DebugData data7({7, 7, 7}, 15);
        DebugData data8({8, 8, 8}, 5);
        auto result = evaluate(data, data1, data2, data3, data4, data5, data6, data7, data8);
        CHECK(result.distances.first   == std::round(1*width));
        CHECK(result.distances.second  == std::round(std::sqrt(3)*width));
        CHECK(result.distances.third   == std::round(std::sqrt(12)*width));
        CHECK(result.distances.fourth  == std::round(std::sqrt(27)*width));
        CHECK(result.distances.fifth   == std::round(std::sqrt(48)*width));
        CHECK(result.distances.sixth   == std::round(std::sqrt(75)*width));
        CHECK(result.distances.seventh == std::round(std::sqrt(108)*width));
        CHECK(result.distances.eighth  == std::round(std::sqrt(147)*width));
        CHECK(result.weights.first == 8);
        CHECK(result.weights.second == 16);
        CHECK(result.weights.third == 32);
        CHECK(result.weights.fourth == 64);
        CHECK(result.weights.fifth == 128);
        CHECK(result.weights.sixth == 256);
        CHECK(result.weights.seventh == 30);
        CHECK(result.weights.eighth == 10);
    }
}

TEST_CASE("CompactCoordinatesData::evaluate") {
    CompactCoordinatesData::inv_width = GENERATE(1, 2, 5);
    SECTION("scalar") {
        // single_tests([](const DebugData& data1, const DebugData& data2) { return data1.evaluate_scalar(data2); });
        // quad_tests([](const DebugData& data, const DebugData& data1, const DebugData& data2, const DebugData& data3, const DebugData& data4) { return data.evaluate_scalar(data1, data2, data3, data4); });
        // octo_tests([](const DebugData& data, const DebugData& data1, const DebugData& data2, const DebugData& data3, const DebugData& data4, const DebugData& data5, const DebugData& data6, const DebugData& data7, const DebugData& data8) { return data.evaluate_scalar(data1, data2, data3, data4, data5, data6, data7, data8); });
    }

    #if defined __SSE2__
        SECTION("sse") {
            // single_tests([](const DebugData& data1, const DebugData& data2) { return data1.evaluate_sse(data2); });
            // quad_tests([](const DebugData& data, const DebugData& data1, const DebugData& data2, const DebugData& data3, const DebugData& data4) { return data.evaluate_sse(data1, data2, data3, data4); });
            // octo_tests([](const DebugData& data, const DebugData& data1, const DebugData& data2, const DebugData& data3, const DebugData& data4, const DebugData& data5, const DebugData& data6, const DebugData& data7, const DebugData& data8) { return data.evaluate_sse(data1, data2, data3, data4, data5, data6, data7, data8); });
        }
    #endif

    #if defined __AVX__
        SECTION("avx") {
            // single_tests([](const DebugData& data1, const DebugData& data2) { return data1.evaluate_avx(data2); });
            quad_tests([](const DebugData& data, const DebugData& data1, const DebugData& data2, const DebugData& data3, const DebugData& data4) { return data.evaluate_avx(data1, data2, data3, data4); });
            octo_tests([](const DebugData& data, const DebugData& data1, const DebugData& data2, const DebugData& data3, const DebugData& data4, const DebugData& data5, const DebugData& data6, const DebugData& data7, const DebugData& data8) { return data.evaluate_avx(data1, data2, data3, data4, data5, data6, data7, data8); });
        }
    #endif
}