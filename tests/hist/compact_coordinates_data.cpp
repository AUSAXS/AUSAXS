#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <hist/detail/CompactCoordinatesData.h>
#include <constants/Constants.h>
#include <math/Vector3.h>

using namespace ausaxs;
using namespace hist::detail;

TEST_CASE("CompactCoordinatesData::CompactCoordinatesData") {
    SECTION("Vector3<double>, float") {
        CompactCoordinatesData data(Vector3<double>(1, 2, 3), 4);
        CHECK(data.value.pos.x() == 1);
        CHECK(data.value.pos.y() == 2);
        CHECK(data.value.pos.z() == 3);
        CHECK(data.value.w == 4);
    }
}

struct DebugData : CompactCoordinatesData {
    using CompactCoordinatesData::CompactCoordinatesData;
    EvaluatedResult evaluate_scalar(const CompactCoordinatesData& other) const {return CompactCoordinatesData::evaluate_scalar(other);}
    QuadEvaluatedResult evaluate_scalar(const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4) const {return CompactCoordinatesData::evaluate_scalar(v1, v2, v3, v4);}
    OctoEvaluatedResult evaluate_scalar(const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4, const CompactCoordinatesData& v5, const CompactCoordinatesData& v6, const CompactCoordinatesData& v7, const CompactCoordinatesData& v8) const {return CompactCoordinatesData::evaluate_scalar(v1, v2, v3, v4, v5, v6, v7, v8);}

    EvaluatedResultRounded evaluate_rounded_scalar(const CompactCoordinatesData& other) const {return CompactCoordinatesData::evaluate_rounded_scalar(other);}
    QuadEvaluatedResultRounded evaluate_rounded_scalar(const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4) const {return CompactCoordinatesData::evaluate_rounded_scalar(v1, v2, v3, v4);}
    OctoEvaluatedResultRounded evaluate_rounded_scalar(const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4, const CompactCoordinatesData& v5, const CompactCoordinatesData& v6, const CompactCoordinatesData& v7, const CompactCoordinatesData& v8) const {return CompactCoordinatesData::evaluate_rounded_scalar(v1, v2, v3, v4, v5, v6, v7, v8);}    

    #if defined __SSE2__
        EvaluatedResult evaluate_sse(const CompactCoordinatesData& other) const {return CompactCoordinatesData::evaluate_sse(other);}
        QuadEvaluatedResult evaluate_sse(const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4) const {return CompactCoordinatesData::evaluate_sse(v1, v2, v3, v4);}
        OctoEvaluatedResult evaluate_sse(const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4, const CompactCoordinatesData& v5, const CompactCoordinatesData& v6, const CompactCoordinatesData& v7, const CompactCoordinatesData& v8) const {return CompactCoordinatesData::evaluate_sse(v1, v2, v3, v4, v5, v6, v7, v8);}

        EvaluatedResultRounded evaluate_rounded_sse(const CompactCoordinatesData& other) const {return CompactCoordinatesData::evaluate_rounded_sse(other);}
        QuadEvaluatedResultRounded evaluate_rounded_sse(const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4) const {return CompactCoordinatesData::evaluate_rounded_sse(v1, v2, v3, v4);}
        OctoEvaluatedResultRounded evaluate_rounded_sse(const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4, const CompactCoordinatesData& v5, const CompactCoordinatesData& v6, const CompactCoordinatesData& v7, const CompactCoordinatesData& v8) const {return CompactCoordinatesData::evaluate_rounded_sse(v1, v2, v3, v4, v5, v6, v7, v8);}
    #endif

    #if defined __AVX__
        EvaluatedResult evaluate_avx(const CompactCoordinatesData& other) const {return CompactCoordinatesData::evaluate_avx(other);}
        QuadEvaluatedResult evaluate_avx(const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4) const {return CompactCoordinatesData::evaluate_avx(v1, v2, v3, v4);}
        OctoEvaluatedResult evaluate_avx(const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4, const CompactCoordinatesData& v5, const CompactCoordinatesData& v6, const CompactCoordinatesData& v7, const CompactCoordinatesData& v8) const {return CompactCoordinatesData::evaluate_avx(v1, v2, v3, v4, v5, v6, v7, v8);}

        EvaluatedResultRounded evaluate_rounded_avx(const CompactCoordinatesData& other) const {return CompactCoordinatesData::evaluate_rounded_avx(other);}
        QuadEvaluatedResultRounded evaluate_rounded_avx(const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4) const {return CompactCoordinatesData::evaluate_rounded_avx(v1, v2, v3, v4);}
        OctoEvaluatedResultRounded evaluate_rounded_avx(const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4, const CompactCoordinatesData& v5, const CompactCoordinatesData& v6, const CompactCoordinatesData& v7, const CompactCoordinatesData& v8) const {return CompactCoordinatesData::evaluate_rounded_avx(v1, v2, v3, v4, v5, v6, v7, v8);}
    #endif
};

#include <functional>
void single_tests(std::function<EvaluatedResult(const DebugData&, const DebugData&)> evaluate) {
    SECTION("single distance") {
        DebugData data1(Vector3<double>{1, 1, 1}, 2);
        DebugData data2(Vector3<double>{2, 1, 1}, 4);
        auto result = evaluate(data1, data2);
        CHECK(result.distance == 1);
        CHECK(result.weight == 8);

        DebugData data3(Vector3<double>{2, 2, 2}, 8);
        result = evaluate(data1, data3);
        CHECK_THAT(result.distance, Catch::Matchers::WithinAbs(std::sqrt(3), 1e-6));
        CHECK(result.weight == 16);
    }
}

void single_tests_rounded(std::function<EvaluatedResultRounded(const DebugData&, const DebugData&)> evaluate) {
    SECTION("single distance") {
        double width = constants::axes::d_axis.width();
        DebugData data1(Vector3<double>{1, 1, 1}, 2);
        DebugData data2(Vector3<double>{2, 1, 1}, 4);
        auto result = evaluate(data1, data2);
        CHECK(result.distance == std::round(1./width));
        CHECK(result.weight == 8);

        DebugData data3(Vector3<double>{2, 2, 2}, 8);
        result = evaluate(data1, data3);
        CHECK(result.distance == std::round(std::sqrt(3)/width));
        CHECK(result.weight == 16);
    }
}

void quad_tests(std::function<QuadEvaluatedResult(const DebugData&, const DebugData&, const DebugData&, const DebugData&, const DebugData&)> evaluate) {
    SECTION("four distances") {
        DebugData data( Vector3<double>{1, 1, 1}, 2);
        DebugData data1(Vector3<double>{2, 1, 1}, 4);
        DebugData data2(Vector3<double>{2, 2, 2}, 8);
        DebugData data3(Vector3<double>{3, 3, 3}, 16);
        DebugData data4(Vector3<double>{4, 4, 4}, 3);
        auto result = evaluate(data, data1, data2, data3, data4);
        CHECK_THAT(result.distances[0],  Catch::Matchers::WithinAbs(1, 1e-6));
        CHECK_THAT(result.distances[1], Catch::Matchers::WithinAbs(std::sqrt(3), 1e-6));
        CHECK_THAT(result.distances[2],  Catch::Matchers::WithinAbs(std::sqrt(12), 1e-6));
        CHECK_THAT(result.distances[3], Catch::Matchers::WithinAbs(std::sqrt(27), 1e-6));
        CHECK(result.weights[0] == 8);
        CHECK(result.weights[1] == 16);
        CHECK(result.weights[2] == 32);
        CHECK(result.weights[3] == 6);
    }
}

void quad_tests_rounded(std::function<QuadEvaluatedResultRounded(const DebugData&, const DebugData&, const DebugData&, const DebugData&, const DebugData&)> evaluate) {
    SECTION("four distances") {
        double width = constants::axes::d_axis.width();
        DebugData data( Vector3<double>{1, 1, 1}, 2);
        DebugData data1(Vector3<double>{2, 1, 1}, 4);
        DebugData data2(Vector3<double>{2, 2, 2}, 8);
        DebugData data3(Vector3<double>{3, 3, 3}, 16);
        DebugData data4(Vector3<double>{4, 4, 4}, 3);
        auto result = evaluate(data, data1, data2, data3, data4);
        CHECK(result.distances[0]  == std::round(1./width));
        CHECK(result.distances[1] == std::round(std::sqrt(3)/width));
        CHECK(result.distances[2]  == std::round(std::sqrt(12)/width));
        CHECK(result.distances[3] == std::round(std::sqrt(27)/width));
        CHECK(result.weights[0] == 8);
        CHECK(result.weights[1] == 16);
        CHECK(result.weights[2] == 32);
        CHECK(result.weights[3] == 6);
    }
}

void octo_tests(std::function<OctoEvaluatedResult(const DebugData&, const DebugData&, const DebugData&, const DebugData&, const DebugData&, const DebugData&, const DebugData&, const DebugData&, const DebugData&)> evaluate) {
    SECTION("eight distances") {
        DebugData data( Vector3<double>{1, 1, 1}, 2);
        DebugData data1(Vector3<double>{2, 1, 1}, 4);
        DebugData data2(Vector3<double>{2, 2, 2}, 8);
        DebugData data3(Vector3<double>{3, 3, 3}, 16);
        DebugData data4(Vector3<double>{4, 4, 4}, 32);
        DebugData data5(Vector3<double>{5, 5, 5}, 64);
        DebugData data6(Vector3<double>{6, 6, 6}, 128);
        DebugData data7(Vector3<double>{7, 7, 7}, 15);
        DebugData data8(Vector3<double>{8, 8, 8}, 5);
        auto result = evaluate(data, data1, data2, data3, data4, data5, data6, data7, data8);
        CHECK_THAT(result.distances[0],  Catch::Matchers::WithinAbs(1, 1e-6));
        CHECK_THAT(result.distances[1], Catch::Matchers::WithinAbs(std::sqrt(3), 1e-6));
        CHECK_THAT(result.distances[2],  Catch::Matchers::WithinAbs(std::sqrt(12), 1e-6));
        CHECK_THAT(result.distances[3], Catch::Matchers::WithinAbs(std::sqrt(27), 1e-6));
        CHECK_THAT(result.distances[4],  Catch::Matchers::WithinAbs(std::sqrt(48), 1e-6));
        CHECK_THAT(result.distances[5],  Catch::Matchers::WithinAbs(std::sqrt(75), 1e-6));
        CHECK_THAT(result.distances[6],Catch::Matchers::WithinAbs(std::sqrt(108), 1e-6));
        CHECK_THAT(result.distances[7], Catch::Matchers::WithinAbs(std::sqrt(147), 1e-5)); // float accuracy loss
        CHECK(result.weights[0] == 8);
        CHECK(result.weights[1] == 16);
        CHECK(result.weights[2] == 32);
        CHECK(result.weights[3] == 64);
        CHECK(result.weights[4] == 128);
        CHECK(result.weights[5] == 256);
        CHECK(result.weights[6] == 30);
        CHECK(result.weights[7] == 10);
    }
}

void octo_tests_rounded(std::function<OctoEvaluatedResultRounded(const DebugData&, const DebugData&, const DebugData&, const DebugData&, const DebugData&, const DebugData&, const DebugData&, const DebugData&, const DebugData&)> evaluate) {
    SECTION("eight distances") {
        double width = constants::axes::d_axis.width();
        DebugData data( Vector3<double>{1, 1, 1}, 2);
        DebugData data1(Vector3<double>{2, 1, 1}, 4);
        DebugData data2(Vector3<double>{2, 2, 2}, 8);
        DebugData data3(Vector3<double>{3, 3, 3}, 16);
        DebugData data4(Vector3<double>{4, 4, 4}, 32);
        DebugData data5(Vector3<double>{5, 5, 5}, 64);
        DebugData data6(Vector3<double>{6, 6, 6}, 128);
        DebugData data7(Vector3<double>{7, 7, 7}, 15);
        DebugData data8(Vector3<double>{8, 8, 8}, 5);
        auto result = evaluate(data, data1, data2, data3, data4, data5, data6, data7, data8);
        CHECK(result.distances[0]   == std::round(1./width));
        CHECK(result.distances[1]  == std::round(std::sqrt(3)/width));
        CHECK(result.distances[2]   == std::round(std::sqrt(12)/width));
        CHECK(result.distances[3]  == std::round(std::sqrt(27)/width));
        CHECK(result.distances[4]   == std::round(std::sqrt(48)/width));
        CHECK(result.distances[5]   == std::round(std::sqrt(75)/width));
        CHECK(result.distances[6] == std::round(std::sqrt(108)/width));
        CHECK(result.distances[7]  == std::round(std::sqrt(147)/width));
        CHECK(result.weights[0] == 8);
        CHECK(result.weights[1] == 16);
        CHECK(result.weights[2] == 32);
        CHECK(result.weights[3] == 64);
        CHECK(result.weights[4] == 128);
        CHECK(result.weights[5] == 256);
        CHECK(result.weights[6] == 30);
        CHECK(result.weights[7] == 10);
    }
}

TEST_CASE("CompactCoordinatesData::evaluate") {
    SECTION("scalar") {
        single_tests([](const DebugData& data1, const DebugData& data2) { return data1.evaluate_scalar(data2); });
        quad_tests([](const DebugData& data, const DebugData& data1, const DebugData& data2, const DebugData& data3, const DebugData& data4) { return data.evaluate_scalar(data1, data2, data3, data4); });
        octo_tests([](const DebugData& data, const DebugData& data1, const DebugData& data2, const DebugData& data3, const DebugData& data4, const DebugData& data5, const DebugData& data6, const DebugData& data7, const DebugData& data8) { return data.evaluate_scalar(data1, data2, data3, data4, data5, data6, data7, data8); });

        single_tests_rounded([](const DebugData& data1, const DebugData& data2) { return data1.evaluate_rounded_scalar(data2); });
        quad_tests_rounded([](const DebugData& data, const DebugData& data1, const DebugData& data2, const DebugData& data3, const DebugData& data4) { return data.evaluate_rounded_scalar(data1, data2, data3, data4); });
        octo_tests_rounded([](const DebugData& data, const DebugData& data1, const DebugData& data2, const DebugData& data3, const DebugData& data4, const DebugData& data5, const DebugData& data6, const DebugData& data7, const DebugData& data8) { return data.evaluate_rounded_scalar(data1, data2, data3, data4, data5, data6, data7, data8); });
    }

    #if defined __SSE2__
        SECTION("sse") {
            single_tests([](const DebugData& data1, const DebugData& data2) { return data1.evaluate_sse(data2); });
            quad_tests([](const DebugData& data, const DebugData& data1, const DebugData& data2, const DebugData& data3, const DebugData& data4) { return data.evaluate_sse(data1, data2, data3, data4); });
            octo_tests([](const DebugData& data, const DebugData& data1, const DebugData& data2, const DebugData& data3, const DebugData& data4, const DebugData& data5, const DebugData& data6, const DebugData& data7, const DebugData& data8) { return data.evaluate_sse(data1, data2, data3, data4, data5, data6, data7, data8); });

            single_tests_rounded([](const DebugData& data1, const DebugData& data2) { return data1.evaluate_rounded_sse(data2); });
            quad_tests_rounded([](const DebugData& data, const DebugData& data1, const DebugData& data2, const DebugData& data3, const DebugData& data4) { return data.evaluate_rounded_sse(data1, data2, data3, data4); });
            octo_tests_rounded([](const DebugData& data, const DebugData& data1, const DebugData& data2, const DebugData& data3, const DebugData& data4, const DebugData& data5, const DebugData& data6, const DebugData& data7, const DebugData& data8) { return data.evaluate_rounded_sse(data1, data2, data3, data4, data5, data6, data7, data8); });
        }
    #endif

    #if defined __AVX__
        SECTION("avx") {
            single_tests([](const DebugData& data1, const DebugData& data2) { return data1.evaluate_avx(data2); });
            quad_tests([](const DebugData& data, const DebugData& data1, const DebugData& data2, const DebugData& data3, const DebugData& data4) { return data.evaluate_avx(data1, data2, data3, data4); });
            octo_tests([](const DebugData& data, const DebugData& data1, const DebugData& data2, const DebugData& data3, const DebugData& data4, const DebugData& data5, const DebugData& data6, const DebugData& data7, const DebugData& data8) { return data.evaluate_avx(data1, data2, data3, data4, data5, data6, data7, data8); });

            single_tests_rounded([](const DebugData& data1, const DebugData& data2) { return data1.evaluate_rounded_avx(data2); });
            quad_tests_rounded([](const DebugData& data, const DebugData& data1, const DebugData& data2, const DebugData& data3, const DebugData& data4) { return data.evaluate_rounded_avx(data1, data2, data3, data4); });
            octo_tests_rounded([](const DebugData& data, const DebugData& data1, const DebugData& data2, const DebugData& data3, const DebugData& data4, const DebugData& data5, const DebugData& data6, const DebugData& data7, const DebugData& data8) { return data.evaluate_rounded_avx(data1, data2, data3, data4, data5, data6, data7, data8); });
        }
    #endif
}