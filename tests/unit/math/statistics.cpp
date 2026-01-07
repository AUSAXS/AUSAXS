#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <math/Statistics.h>

using namespace ausaxs;

TEST_CASE("stats::mean") {
    SECTION("simple cases") {
        CHECK(stats::mean(std::vector{10, 3, 5, 6}) == 6);
        CHECK(stats::mean(std::vector{12, 14, 15, 15, 14, 17}) == 14.5);
        CHECK(stats::mean(std::vector{54, 66, 78, 80, 82, 84, 84, 90, 93}) == 79);
    }

    SECTION("single element") {
        CHECK(stats::mean(std::vector{42}) == 42);
    }

    SECTION("negative numbers") {
        CHECK(stats::mean(std::vector{-10, -5, 0, 5, 10}) == 0);
    }

    SECTION("floating point") {
        CHECK_THAT(stats::mean(std::vector{1.5, 2.5, 3.5}), Catch::Matchers::WithinAbs(2.5, 1e-10));
    }
}

TEST_CASE("stats::var") {
    SECTION("simple cases") {
        CHECK(stats::var(std::vector{9, 10, 11, 7, 13}) == 5);
        CHECK(stats::var(std::vector{10, 10, 10, 10, 10}) == 0);
        CHECK(stats::var(std::vector{1, 1, 10, 19, 19}) == 81);
    }

    SECTION("ddof parameter") {
        std::vector<double> data = {1, 2, 3, 4, 5};
        double var1 = stats::var(data, 1);
        double var0 = stats::var(data, 0);
        CHECK(var0 < var1);
    }
}

TEST_CASE("stats::std") {
    SECTION("simple cases") {
        CHECK_THAT(stats::std(std::vector{9, 10, 11, 7, 13}), Catch::Matchers::WithinRel(std::sqrt(5)));
        CHECK(stats::std(std::vector{10, 10, 10, 10, 10}) == 0);
        CHECK(stats::std(std::vector{1, 1, 10, 19, 19}) == 9);
    }

    SECTION("relationship to variance") {
        std::vector<double> data = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
        double variance = stats::var(data);
        double std_dev = stats::std(data);
        CHECK_THAT(std_dev * std_dev, Catch::Matchers::WithinAbs(variance, 1e-10));
    }
}

TEST_CASE("stats::mode") {
    SECTION("simple cases") {
        CHECK(stats::mode(std::vector{10, 3, 5, 6, 6}) == 6);
        CHECK(stats::mode(std::vector{12, 14, 15, 15, 14, 17, 15}) == 15);
        CHECK(stats::mode(std::vector{54, 66, 78, 80, 82, 84, 84, 90, 93}) == 84);
        CHECK(stats::mode(std::vector{8, 1, 3, 5, 8, 6, 3, 7, 8, 9, 3, 2, 1}) == 8);
    }

    SECTION("single element") {
        CHECK(stats::mode(std::vector{42}) == 42);
    }

    SECTION("all unique") {
        int result = stats::mode(std::vector{1, 2, 3, 4, 5});
        CHECK((result >= 1 && result <= 5));
    }
}

TEST_CASE("stats::weighted_mean") {
    SECTION("uniform weights") {
        auto v = std::vector{10, 3, 5, 6};
        CHECK_THAT(stats::weighted_mean(v, std::vector{1, 1, 1, 1}), Catch::Matchers::WithinAbs(stats::mean(v), 1e-3));
    }

    SECTION("non-uniform weights") {
        auto v = std::vector{10, 3, 5, 6};
        CHECK_THAT(stats::weighted_mean(v, std::vector{1, 2, 3, 4}), Catch::Matchers::WithinAbs(8.204903, 1e-3));
    }

    SECTION("more complex case") {
        auto v = std::vector{12, 14, 15, 15, 14, 17};
        CHECK_THAT(stats::weighted_mean(v, std::vector<double>{1, 1, 1, 1, 1, 1}), Catch::Matchers::WithinAbs(stats::mean(v), 1e-3));
        CHECK_THAT(stats::weighted_mean(v, std::vector<double>{1, 0.5, 0.1, 2, 0.9, 3}), Catch::Matchers::WithinAbs(14.924834, 1e-3));
    }

    SECTION("large dataset") {
        auto v = std::vector{54, 66, 78, 80, 82, 84, 84, 90, 93};
        CHECK_THAT(stats::weighted_mean(v, std::vector<double>{1, 1, 1, 1, 1, 1, 1, 1, 1}), Catch::Matchers::WithinAbs(stats::mean(v), 1e-3));       
        CHECK_THAT(stats::weighted_mean(v, std::vector<double>{1, 0.2, 2, 0.5, 0.9, 4, 1, 0.1, 0.4}), Catch::Matchers::WithinAbs(85.125966, 1e-3));       
    }
}

TEST_CASE("stats::weighted_mean_error") {
    SECTION("simple cases") {
        CHECK_THAT(stats::weighted_mean_error(std::vector{1, 2, 3, 4}), Catch::Matchers::WithinAbs(0.838116, 1e-6));
        CHECK_THAT(stats::weighted_mean_error(std::vector<double>{1, 0.5, 0.1, 2, 0.9, 3}), Catch::Matchers::WithinAbs(0.0968561, 1e-6));
        CHECK_THAT(stats::weighted_mean_error(std::vector<double>{1, 0.2, 2, 0.5, 0.9, 4, 1, 0.1, 0.4}), Catch::Matchers::WithinAbs(0.0848808, 1e-6));
    }
}

TEST_CASE("stats::Measurement") {
    SECTION("construction") {
        std::vector<double> data = {1, 2, 3, 4, 5};
        stats::Measurement<double> m(data);
        CHECK(m.vals == data);
    }

    SECTION("mean") {
        std::vector<double> data = {10, 3, 5, 6};
        stats::Measurement<double> m(data);
        CHECK(stats::mean(m.vals) == 6);
    }

    SECTION("variance") {
        std::vector<double> data = {9, 10, 11, 7, 13};
        stats::Measurement<double> m(data);
        CHECK(stats::var(m.vals) == 5);
    }

    SECTION("standard deviation") {
        std::vector<double> data = {10, 10, 10, 10, 10};
        stats::Measurement<double> m(data);
        CHECK(stats::std(m.vals) == 0);
    }
}
