#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <math/MovingAverager.h>

using namespace ausaxs;

TEST_CASE("MovingAverage::average") {
    SECTION("simple uniform data") {
        std::vector<double> data = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
        
        SECTION("window size 3") {
            auto res = MovingAverage::average(data, 3);
            REQUIRE(res.size() == 10);
            CHECK(res[0] == 1);
            CHECK(res[1] == 2);
            CHECK(res[2] == 3);
            CHECK(res[3] == 4);
            CHECK(res[4] == 5);
            CHECK(res[5] == 6);
            CHECK(res[6] == 7);
            CHECK(res[7] == 8);
            CHECK(res[8] == 9);
            CHECK(res[9] == 10);
        }

        SECTION("window size 5") {
            auto res = MovingAverage::average(data, 5);
            REQUIRE(res.size() == 10);
            CHECK(res[0] == 1);
            CHECK(res[1] == 2);
            CHECK(res[2] == 3);
            CHECK(res[3] == 4);
            CHECK(res[4] == 5);
            CHECK(res[5] == 6);
            CHECK(res[6] == 7);
            CHECK(res[7] == 8);
            CHECK(res[8] == 9);
            CHECK(res[9] == 10);
        }
    }

    SECTION("varying data") {
        std::vector<double> data = {1, 5, 3, 7, 2, 8, 4, 6};
        
        SECTION("window size 3") {
            auto res = MovingAverage::average(data, 3);
            REQUIRE(res.size() == 8);
            CHECK_THAT(res[1], Catch::Matchers::WithinAbs(3.0, 1e-10));
        }
    }

    SECTION("edge cases preserved") {
        std::vector<double> data = {10, 20, 30, 40, 50};
        auto res = MovingAverage::average(data, 3);
        REQUIRE(res.size() == 5);
        CHECK_THAT(res[0], Catch::Matchers::WithinAbs(10.0, 1e-10));
        CHECK_THAT(res[4], Catch::Matchers::WithinAbs(50.0, 1e-10));
    }
}

TEST_CASE("MovingAverage::average_half") {
    SECTION("simple data") {
        std::vector<double> data = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
        
        auto res = MovingAverage::average_half(data, 3);
        REQUIRE(res.size() == 10);
    }

    SECTION("weighted averaging") {
        std::vector<double> data = {10, 20, 30, 40, 50};
        auto res = MovingAverage::average_half(data, 3);
        REQUIRE(res.size() == 5);
    }
}
