#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>

#include <math/Statistics.h>

using std::vector;

TEST_CASE("basic_stats", "[stats]") {
    SECTION("mean") {
        CHECK(stats::mean({10, 3, 5, 6}) == 6);
        CHECK(stats::mean({12, 14, 15, 15, 14, 17}) == 14.5);
        CHECK(stats::mean({54, 66, 78, 80, 82, 84, 84, 90, 93}) == 79);
    }

    SECTION("std") {
        CHECK_THAT(stats::std({9, 10, 11, 7, 13}), Catch::Matchers::WithinRel(std::sqrt(5)));
        CHECK(stats::std({10, 10, 10, 10, 10}) == 0);
        CHECK(stats::std({1, 1, 10, 19, 19}) == 9);
    }

    SECTION("var") {
        CHECK(stats::var({9, 10, 11, 7, 13}) == 5);
        CHECK(stats::var({10, 10, 10, 10, 10}) == 0);
        CHECK(stats::var({1, 1, 10, 19, 19}) == 81);
    }
}

TEST_CASE("math_mode", "[stats]") {
    SECTION("mode") {
        CHECK(stats::mode(vector{10, 3, 5, 6, 6}) == 6);
        CHECK(stats::mode(vector{12, 14, 15, 15, 14, 17, 15}) == 15);
        CHECK(stats::mode(vector{54, 66, 78, 80, 82, 84, 84, 90, 93}) == 84);
        CHECK(stats::mode(vector{8, 1, 3, 5, 8, 6, 3, 7, 8, 9, 3, 2, 1}) == 8);
    }
}