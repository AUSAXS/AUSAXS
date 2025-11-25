#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <hist/detail/data/CompactCoordinatesXYZW.h>
#include <constants/Constants.h>
#include <math/Vector3.h>

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
    using CompactCoordinatesXYZW<vbw>::CompactCoordinatesXYZW;
    EvaluatedResult evaluate_scalar(const CompactCoordinatesXYZW<vbw>& other) const {return CompactCoordinatesXYZW<vbw>::evaluate_scalar(other);}
    QuadEvaluatedResult evaluate_scalar(const CompactCoordinatesXYZW<vbw>& v1, const CompactCoordinatesXYZW<vbw>& v2, const CompactCoordinatesXYZW<vbw>& v3, const CompactCoordinatesXYZW<vbw>& v4) const {return CompactCoordinatesXYZW<vbw>::evaluate_scalar(v1, v2, v3, v4);}
    OctoEvaluatedResult evaluate_scalar(const CompactCoordinatesXYZW<vbw>& v1, const CompactCoordinatesXYZW<vbw>& v2, const CompactCoordinatesXYZW<vbw>& v3, const CompactCoordinatesXYZW<vbw>& v4, const CompactCoordinatesXYZW<vbw>& v5, const CompactCoordinatesXYZW<vbw>& v6, const CompactCoordinatesXYZW<vbw>& v7, const CompactCoordinatesXYZW<vbw>& v8) const {return CompactCoordinatesXYZW<vbw>::evaluate_scalar(v1, v2, v3, v4, v5, v6, v7, v8);}

    EvaluatedResultRounded evaluate_rounded_scalar(const CompactCoordinatesXYZW<vbw>& other) const {return CompactCoordinatesXYZW<vbw>::evaluate_rounded_scalar(other);}
    QuadEvaluatedResultRounded evaluate_rounded_scalar(const CompactCoordinatesXYZW<vbw>& v1, const CompactCoordinatesXYZW<vbw>& v2, const CompactCoordinatesXYZW<vbw>& v3, const CompactCoordinatesXYZW<vbw>& v4) const {return CompactCoordinatesXYZW<vbw>::evaluate_rounded_scalar(v1, v2, v3, v4);}
    OctoEvaluatedResultRounded evaluate_rounded_scalar(const CompactCoordinatesXYZW<vbw>& v1, const CompactCoordinatesXYZW<vbw>& v2, const CompactCoordinatesXYZW<vbw>& v3, const CompactCoordinatesXYZW<vbw>& v4, const CompactCoordinatesXYZW<vbw>& v5, const CompactCoordinatesXYZW<vbw>& v6, const CompactCoordinatesXYZW<vbw>& v7, const CompactCoordinatesXYZW<vbw>& v8) const {return CompactCoordinatesXYZW<vbw>::evaluate_rounded_scalar(v1, v2, v3, v4, v5, v6, v7, v8);}    

    #if defined __SSE2__
        EvaluatedResult evaluate_sse(const CompactCoordinatesXYZW<vbw>& other) const {return CompactCoordinatesXYZW<vbw>::evaluate_sse(other);}
        QuadEvaluatedResult evaluate_sse(const CompactCoordinatesXYZW<vbw>& v1, const CompactCoordinatesXYZW<vbw>& v2, const CompactCoordinatesXYZW<vbw>& v3, const CompactCoordinatesXYZW<vbw>& v4) const {return CompactCoordinatesXYZW<vbw>::evaluate_sse(v1, v2, v3, v4);}
        OctoEvaluatedResult evaluate_sse(const CompactCoordinatesXYZW<vbw>& v1, const CompactCoordinatesXYZW<vbw>& v2, const CompactCoordinatesXYZW<vbw>& v3, const CompactCoordinatesXYZW<vbw>& v4, const CompactCoordinatesXYZW<vbw>& v5, const CompactCoordinatesXYZW<vbw>& v6, const CompactCoordinatesXYZW<vbw>& v7, const CompactCoordinatesXYZW<vbw>& v8) const {return CompactCoordinatesXYZW<vbw>::evaluate_sse(v1, v2, v3, v4, v5, v6, v7, v8);}

        EvaluatedResultRounded evaluate_rounded_sse(const CompactCoordinatesXYZW<vbw>& other) const {return CompactCoordinatesXYZW<vbw>::evaluate_rounded_sse(other);}
        QuadEvaluatedResultRounded evaluate_rounded_sse(const CompactCoordinatesXYZW<vbw>& v1, const CompactCoordinatesXYZW<vbw>& v2, const CompactCoordinatesXYZW<vbw>& v3, const CompactCoordinatesXYZW<vbw>& v4) const {return CompactCoordinatesXYZW<vbw>::evaluate_rounded_sse(v1, v2, v3, v4);}
        OctoEvaluatedResultRounded evaluate_rounded_sse(const CompactCoordinatesXYZW<vbw>& v1, const CompactCoordinatesXYZW<vbw>& v2, const CompactCoordinatesXYZW<vbw>& v3, const CompactCoordinatesXYZW<vbw>& v4, const CompactCoordinatesXYZW<vbw>& v5, const CompactCoordinatesXYZW<vbw>& v6, const CompactCoordinatesXYZW<vbw>& v7, const CompactCoordinatesXYZW<vbw>& v8) const {return CompactCoordinatesXYZW<vbw>::evaluate_rounded_sse(v1, v2, v3, v4, v5, v6, v7, v8);}
    #endif

    #if defined __AVX__
        EvaluatedResult evaluate_avx(const CompactCoordinatesXYZW<vbw>& other) const {return CompactCoordinatesXYZW<vbw>::evaluate_avx(other);}
        QuadEvaluatedResult evaluate_avx(const CompactCoordinatesXYZW<vbw>& v1, const CompactCoordinatesXYZW<vbw>& v2, const CompactCoordinatesXYZW<vbw>& v3, const CompactCoordinatesXYZW<vbw>& v4) const {return CompactCoordinatesXYZW<vbw>::evaluate_avx(v1, v2, v3, v4);}
        OctoEvaluatedResult evaluate_avx(const CompactCoordinatesXYZW<vbw>& v1, const CompactCoordinatesXYZW<vbw>& v2, const CompactCoordinatesXYZW<vbw>& v3, const CompactCoordinatesXYZW<vbw>& v4, const CompactCoordinatesXYZW<vbw>& v5, const CompactCoordinatesXYZW<vbw>& v6, const CompactCoordinatesXYZW<vbw>& v7, const CompactCoordinatesXYZW<vbw>& v8) const {return CompactCoordinatesXYZW<vbw>::evaluate_avx(v1, v2, v3, v4, v5, v6, v7, v8);}

        EvaluatedResultRounded evaluate_rounded_avx(const CompactCoordinatesXYZW<vbw>& other) const {return CompactCoordinatesXYZW<vbw>::evaluate_rounded_avx(other);}
        QuadEvaluatedResultRounded evaluate_rounded_avx(const CompactCoordinatesXYZW<vbw>& v1, const CompactCoordinatesXYZW<vbw>& v2, const CompactCoordinatesXYZW<vbw>& v3, const CompactCoordinatesXYZW<vbw>& v4) const {return CompactCoordinatesXYZW<vbw>::evaluate_rounded_avx(v1, v2, v3, v4);}
        OctoEvaluatedResultRounded evaluate_rounded_avx(const CompactCoordinatesXYZW<vbw>& v1, const CompactCoordinatesXYZW<vbw>& v2, const CompactCoordinatesXYZW<vbw>& v3, const CompactCoordinatesXYZW<vbw>& v4, const CompactCoordinatesXYZW<vbw>& v5, const CompactCoordinatesXYZW<vbw>& v6, const CompactCoordinatesXYZW<vbw>& v7, const CompactCoordinatesXYZW<vbw>& v8) const {return CompactCoordinatesXYZW<vbw>::evaluate_rounded_avx(v1, v2, v3, v4, v5, v6, v7, v8);}
    #endif
};

#include <functional>
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
void quad_tests(std::function<QuadEvaluatedResult(const DebugData<vbw>&, const DebugData<vbw>&, const DebugData<vbw>&, const DebugData<vbw>&, const DebugData<vbw>&)> evaluate) {
    SECTION("four distances") {
        DebugData<vbw> data( Vector3<double>{1, 1, 1}, 2);
        DebugData<vbw> data1(Vector3<double>{2, 1, 1}, 4);
        DebugData<vbw> data2(Vector3<double>{2, 2, 2}, 8);
        DebugData<vbw> data3(Vector3<double>{3, 3, 3}, 16);
        DebugData<vbw> data4(Vector3<double>{4, 4, 4}, 3);
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

template<bool vbw>
void quad_tests_rounded(std::function<QuadEvaluatedResultRounded(const DebugData<vbw>&, const DebugData<vbw>&, const DebugData<vbw>&, const DebugData<vbw>&, const DebugData<vbw>&)> evaluate) {
    SECTION("four distances") {
        double width = constants::axes::d_axis.width();
        DebugData<vbw> data( Vector3<double>{1, 1, 1}, 2);
        DebugData<vbw> data1(Vector3<double>{2, 1, 1}, 4);
        DebugData<vbw> data2(Vector3<double>{2, 2, 2}, 8);
        DebugData<vbw> data3(Vector3<double>{3, 3, 3}, 16);
        DebugData<vbw> data4(Vector3<double>{4, 4, 4}, 3);
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

template<bool vbw>
void octo_tests(std::function<OctoEvaluatedResult(const DebugData<vbw>&, const DebugData<vbw>&, const DebugData<vbw>&, const DebugData<vbw>&, const DebugData<vbw>&, const DebugData<vbw>&, const DebugData<vbw>&, const DebugData<vbw>&, const DebugData<vbw>&)> evaluate) {
    SECTION("eight distances") {
        DebugData<vbw> data( Vector3<double>{1, 1, 1}, 2);
        DebugData<vbw> data1(Vector3<double>{2, 1, 1}, 4);
        DebugData<vbw> data2(Vector3<double>{2, 2, 2}, 8);
        DebugData<vbw> data3(Vector3<double>{3, 3, 3}, 16);
        DebugData<vbw> data4(Vector3<double>{4, 4, 4}, 32);
        DebugData<vbw> data5(Vector3<double>{5, 5, 5}, 64);
        DebugData<vbw> data6(Vector3<double>{6, 6, 6}, 128);
        DebugData<vbw> data7(Vector3<double>{7, 7, 7}, 15);
        DebugData<vbw> data8(Vector3<double>{8, 8, 8}, 5);
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

template<bool vbw>
void octo_tests_rounded(std::function<OctoEvaluatedResultRounded(const DebugData<vbw>&, const DebugData<vbw>&, const DebugData<vbw>&, const DebugData<vbw>&, const DebugData<vbw>&, const DebugData<vbw>&, const DebugData<vbw>&, const DebugData<vbw>&, const DebugData<vbw>&)> evaluate) {
    SECTION("eight distances") {
        double width = constants::axes::d_axis.width();
        DebugData<vbw> data( Vector3<double>{1, 1, 1}, 2);
        DebugData<vbw> data1(Vector3<double>{2, 1, 1}, 4);
        DebugData<vbw> data2(Vector3<double>{2, 2, 2}, 8);
        DebugData<vbw> data3(Vector3<double>{3, 3, 3}, 16);
        DebugData<vbw> data4(Vector3<double>{4, 4, 4}, 32);
        DebugData<vbw> data5(Vector3<double>{5, 5, 5}, 64);
        DebugData<vbw> data6(Vector3<double>{6, 6, 6}, 128);
        DebugData<vbw> data7(Vector3<double>{7, 7, 7}, 15);
        DebugData<vbw> data8(Vector3<double>{8, 8, 8}, 5);
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

template<bool vbw>
void run_tests() {
    SECTION("scalar") {
        single_tests<vbw>([](const DebugData<vbw>& data1, const DebugData<vbw>& data2) { return data1.evaluate_scalar(data2); });
        quad_tests<vbw>([](const DebugData<vbw>& data, const DebugData<vbw>& data1, const DebugData<vbw>& data2, const DebugData<vbw>& data3, const DebugData<vbw>& data4) { return data.evaluate_scalar(data1, data2, data3, data4); });
        octo_tests<vbw>([](const DebugData<vbw>& data, const DebugData<vbw>& data1, const DebugData<vbw>& data2, const DebugData<vbw>& data3, const DebugData<vbw>& data4, const DebugData<vbw>& data5, const DebugData<vbw>& data6, const DebugData<vbw>& data7, const DebugData<vbw>& data8) { return data.evaluate_scalar(data1, data2, data3, data4, data5, data6, data7, data8); });

        single_tests_rounded<vbw>([](const DebugData<vbw>& data1, const DebugData<vbw>& data2) { return data1.evaluate_rounded_scalar(data2); });
        quad_tests_rounded<vbw>([](const DebugData<vbw>& data, const DebugData<vbw>& data1, const DebugData<vbw>& data2, const DebugData<vbw>& data3, const DebugData<vbw>& data4) { return data.evaluate_rounded_scalar(data1, data2, data3, data4); });
        octo_tests_rounded<vbw>([](const DebugData<vbw>& data, const DebugData<vbw>& data1, const DebugData<vbw>& data2, const DebugData<vbw>& data3, const DebugData<vbw>& data4, const DebugData<vbw>& data5, const DebugData<vbw>& data6, const DebugData<vbw>& data7, const DebugData<vbw>& data8) { return data.evaluate_rounded_scalar(data1, data2, data3, data4, data5, data6, data7, data8); });
    }

    #if defined __SSE2__
        SECTION("sse") {
            single_tests<vbw>([](const DebugData<vbw>& data1, const DebugData<vbw>& data2) { return data1.evaluate_sse(data2); });
            quad_tests<vbw>([](const DebugData<vbw>& data, const DebugData<vbw>& data1, const DebugData<vbw>& data2, const DebugData<vbw>& data3, const DebugData<vbw>& data4) { return data.evaluate_sse(data1, data2, data3, data4); });
            octo_tests<vbw>([](const DebugData<vbw>& data, const DebugData<vbw>& data1, const DebugData<vbw>& data2, const DebugData<vbw>& data3, const DebugData<vbw>& data4, const DebugData<vbw>& data5, const DebugData<vbw>& data6, const DebugData<vbw>& data7, const DebugData<vbw>& data8) { return data.evaluate_sse(data1, data2, data3, data4, data5, data6, data7, data8); });

            single_tests_rounded<vbw>([](const DebugData<vbw>& data1, const DebugData<vbw>& data2) { return data1.evaluate_rounded_sse(data2); });
            quad_tests_rounded<vbw>([](const DebugData<vbw>& data, const DebugData<vbw>& data1, const DebugData<vbw>& data2, const DebugData<vbw>& data3, const DebugData<vbw>& data4) { return data.evaluate_rounded_sse(data1, data2, data3, data4); });
            octo_tests_rounded<vbw>([](const DebugData<vbw>& data, const DebugData<vbw>& data1, const DebugData<vbw>& data2, const DebugData<vbw>& data3, const DebugData<vbw>& data4, const DebugData<vbw>& data5, const DebugData<vbw>& data6, const DebugData<vbw>& data7, const DebugData<vbw>& data8) { return data.evaluate_rounded_sse(data1, data2, data3, data4, data5, data6, data7, data8); });
        }
    #endif

    #if defined __AVX__
        SECTION("avx") {
            single_tests<vbw>([](const DebugData<vbw>& data1, const DebugData<vbw>& data2) { return data1.evaluate_avx(data2); });
            quad_tests<vbw>([](const DebugData<vbw>& data, const DebugData<vbw>& data1, const DebugData<vbw>& data2, const DebugData<vbw>& data3, const DebugData<vbw>& data4) { return data.evaluate_avx(data1, data2, data3, data4); });
            octo_tests<vbw>([](const DebugData<vbw>& data, const DebugData<vbw>& data1, const DebugData<vbw>& data2, const DebugData<vbw>& data3, const DebugData<vbw>& data4, const DebugData<vbw>& data5, const DebugData<vbw>& data6, const DebugData<vbw>& data7, const DebugData<vbw>& data8) { return data.evaluate_avx(data1, data2, data3, data4, data5, data6, data7, data8); });

            single_tests_rounded<vbw>([](const DebugData<vbw>& data1, const DebugData<vbw>& data2) { return data1.evaluate_rounded_avx(data2); });
            quad_tests_rounded<vbw>([](const DebugData<vbw>& data, const DebugData<vbw>& data1, const DebugData<vbw>& data2, const DebugData<vbw>& data3, const DebugData<vbw>& data4) { return data.evaluate_rounded_avx(data1, data2, data3, data4); });
            octo_tests_rounded<vbw>([](const DebugData<vbw>& data, const DebugData<vbw>& data1, const DebugData<vbw>& data2, const DebugData<vbw>& data3, const DebugData<vbw>& data4, const DebugData<vbw>& data5, const DebugData<vbw>& data6, const DebugData<vbw>& data7, const DebugData<vbw>& data8) { return data.evaluate_rounded_avx(data1, data2, data3, data4, data5, data6, data7, data8); });
        }
    #endif
}

TEST_CASE("CompactCoordinatesXYZW<vbw>::evaluate") {
    SECTION("variable bin width") {
        run_tests<true>();
    }
    SECTION("fixed bin width") {
        run_tests<false>();
    }
}