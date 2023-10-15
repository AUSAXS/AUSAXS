#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <hist/detail/CompactCoordinatesData.h>
#include <math/Vector3.h>

using namespace hist::detail;

TEST_CASE("CompactCoordinatesData::CompactCoordinatesData") {
    SECTION("default") {
        CompactCoordinatesData data;
        CHECK(data.x == 0);
        CHECK(data.y == 0);
        CHECK(data.z == 0);
        CHECK(data.w == 0);
    }

    SECTION("Vector3<double>, float") {
        CompactCoordinatesData data(Vector3<double>(1, 2, 3), 4);
        CHECK(data.x == 1);
        CHECK(data.y == 2);
        CHECK(data.z == 3);
        CHECK(data.w == 4);
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
    SECTION("single distance") {
        DebugData data1({1, 1, 1}, 2);
        DebugData data2({2, 1, 1}, 4);
        auto result = evaluate(data1, data2);
        CHECK(result.distance == 1);
        CHECK(result.weight == 8);

        DebugData data3({2, 2, 2}, 8);
        result = evaluate(data1, data3);
        CHECK_THAT(result.distance, Catch::Matchers::WithinAbs(std::sqrt(3), 1e-6));
        CHECK(result.weight == 16);
    }
}

void quad_tests(std::function<QuadEvaluatedResult(const DebugData&, const DebugData&, const DebugData&, const DebugData&, const DebugData&)> evaluate) {
    SECTION("four distances") {
        DebugData data({1, 1, 1}, 2);
        DebugData data1({2, 1, 1}, 4);
        DebugData data2({2, 2, 2}, 8);
        DebugData data3({3, 3, 3}, 16);
        DebugData data4({4, 4, 4}, 3);
        auto result = evaluate(data, data1, data2, data3, data4);
        CHECK(     result.distances.first == 1);
        CHECK_THAT(result.distances.second, Catch::Matchers::WithinAbs(std::sqrt(3), 1e-6));
        CHECK_THAT(result.distances.third, Catch::Matchers::WithinAbs(std::sqrt(12), 1e-6));
        CHECK_THAT(result.distances.fourth, Catch::Matchers::WithinAbs(std::sqrt(27), 1e-6));
        CHECK(result.weights.first == 8);
        CHECK(result.weights.second == 16);
        CHECK(result.weights.third == 32);
        CHECK(result.weights.fourth == 6);
    }
}

void octo_tests(std::function<OctoEvaluatedResult(const DebugData&, const DebugData&, const DebugData&, const DebugData&, const DebugData&, const DebugData&, const DebugData&, const DebugData&, const DebugData&)> evaluate) {
    SECTION("eight distances") {
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
        CHECK(     result.distance[0] == 1);
        CHECK_THAT(result.distance[1], Catch::Matchers::WithinAbs(std::sqrt(3), 1e-6));
        CHECK_THAT(result.distance[2], Catch::Matchers::WithinAbs(std::sqrt(12), 1e-6));
        CHECK_THAT(result.distance[3], Catch::Matchers::WithinAbs(std::sqrt(27), 1e-6));
        CHECK_THAT(result.distance[4], Catch::Matchers::WithinAbs(std::sqrt(48), 1e-6));
        CHECK_THAT(result.distance[5], Catch::Matchers::WithinAbs(std::sqrt(75), 1e-6));
        CHECK_THAT(result.distance[6], Catch::Matchers::WithinAbs(std::sqrt(108), 1e-6));
        CHECK_THAT(result.distance[7], Catch::Matchers::WithinAbs(std::sqrt(147), 1e-6));
        CHECK(result.weight[0] == 8);
        CHECK(result.weight[1] == 16);
        CHECK(result.weight[2] == 32);
        CHECK(result.weight[3] == 64);
        CHECK(result.weight[4] == 128);
        CHECK(result.weight[5] == 256);
        CHECK(result.weight[6] == 30);
        CHECK(result.weight[7] == 10);
    }
}

TEST_CASE("CompactCoordinatesData::evaluate") {
    SECTION("scalar") {
        single_tests([](const DebugData& data1, const DebugData& data2) { return data1.evaluate_scalar(data2); });
        quad_tests([](const DebugData& data, const DebugData& data1, const DebugData& data2, const DebugData& data3, const DebugData& data4) { return data.evaluate_scalar(data1, data2, data3, data4); });
        octo_tests([](const DebugData& data, const DebugData& data1, const DebugData& data2, const DebugData& data3, const DebugData& data4, const DebugData& data5, const DebugData& data6, const DebugData& data7, const DebugData& data8) { return data.evaluate_scalar(data1, data2, data3, data4, data5, data6, data7, data8); });
    }

    #if defined __SSE2__
        SECTION("sse") {
            std::cout << "sse enabled" << std::endl;
            single_tests([](const DebugData& data1, const DebugData& data2) { return data1.evaluate_sse(data2); });
            quad_tests([](const DebugData& data, const DebugData& data1, const DebugData& data2, const DebugData& data3, const DebugData& data4) { return data.evaluate_sse(data1, data2, data3, data4); });
            octo_tests([](const DebugData& data, const DebugData& data1, const DebugData& data2, const DebugData& data3, const DebugData& data4, const DebugData& data5, const DebugData& data6, const DebugData& data7, const DebugData& data8) { return data.evaluate_sse(data1, data2, data3, data4, data5, data6, data7, data8); });
        }
    #endif

    #if defined __AVX__
        SECTION("avx") {
            std::cout << "avx enabled" << std::endl;
            // single_tests([](const DebugData& data1, const DebugData& data2) { return data1.evaluate_avx(data2); });
        }
    #endif
}