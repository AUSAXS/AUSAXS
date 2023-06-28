#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>

#include <math/Statistics.h>

using std::vector;

TEST_CASE("basics") {
    SECTION("weighted_mean_error") {
        CHECK_THAT(stats::weighted_mean_error(vector{1, 2, 3, 4}), Catch::Matchers::WithinAbs(0.838116, 1e-6));
        CHECK_THAT(stats::weighted_mean_error(vector<double>{1, 0.5, 0.1, 2, 0.9, 3}), Catch::Matchers::WithinAbs(0.0968561, 1e-6));
        CHECK_THAT(stats::weighted_mean_error(vector<double>{1, 0.2, 2, 0.5, 0.9, 4, 1, 0.1, 0.4}), Catch::Matchers::WithinAbs(0.0848808, 1e-6));
    }

    SECTION("weighted_mean") {
        auto v = vector{10, 3, 5, 6};
        CHECK_THAT(stats::weighted_mean(v, vector{1, 1, 1, 1}), Catch::Matchers::WithinAbs(stats::mean(v), 1e-3));
        CHECK_THAT(stats::weighted_mean(v, vector{1, 2, 3, 4}), Catch::Matchers::WithinAbs(8.204903, 1e-3));

        v = vector{12, 14, 15, 15, 14, 17};
        CHECK_THAT(stats::weighted_mean(v, vector<double>{1, 1, 1, 1, 1, 1}), Catch::Matchers::WithinAbs(stats::mean(v), 1e-3));
        CHECK_THAT(stats::weighted_mean(v, vector<double>{1, 0.5, 0.1, 2, 0.9, 3}), Catch::Matchers::WithinAbs(14.924834, 1e-3));

        v = vector{54, 66, 78, 80, 82, 84, 84, 90, 93};
        CHECK_THAT(stats::weighted_mean(v, vector<double>{1, 1, 1, 1, 1, 1, 1, 1, 1}), Catch::Matchers::WithinAbs(stats::mean(v), 1e-3));       
        CHECK_THAT(stats::weighted_mean(v, vector<double>{1, 0.2, 2, 0.5, 0.9, 4, 1, 0.1, 0.4}), Catch::Matchers::WithinAbs(85.125966, 1e-3));       
    }

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

TEST_CASE("mode") {
    CHECK(stats::mode(vector{10, 3, 5, 6, 6}) == 6);
    CHECK(stats::mode(vector{12, 14, 15, 15, 14, 17, 15}) == 15);
    CHECK(stats::mode(vector{54, 66, 78, 80, 82, 84, 84, 90, 93}) == 84);
    CHECK(stats::mode(vector{8, 1, 3, 5, 8, 6, 3, 7, 8, 9, 3, 2, 1}) == 8);
}