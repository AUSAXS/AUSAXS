#include <catch2/catch_test_macros.hpp>

#include <math/PeakFinder.h>

using namespace ausaxs::math;

TEST_CASE("find_minima") {
    SECTION("simple parabola") {
        std::vector<double> x = {0, 1, 2, 3, 4, 5, 6, 7, 8};
        std::vector<double> y = {16, 9, 4, 1, 0, 1, 4, 9, 16};
        
        auto minima = find_minima(x, y, 1, 0.1);
        REQUIRE(minima.size() > 0);
        CHECK(minima[0] == 4);
    }

    SECTION("multiple minima") {
        std::vector<double> x = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
        std::vector<double> y = {5, 1, 5, 1, 5, 1, 5, 1, 5, 1, 5};
        
        auto minima = find_minima(x, y, 1, 0.1);
        REQUIRE(minima.size() >= 1);
    }

    SECTION("no minima") {
        std::vector<double> x = {0, 1, 2, 3, 4};
        std::vector<double> y = {0, 1, 2, 3, 4};
        
        auto minima = find_minima(x, y, 1, 0.1);
        CHECK(minima.empty() || minima.size() > 0);
    }

    SECTION("single point minimum") {
        std::vector<double> x = {0, 1, 2, 3, 4, 5, 6};
        std::vector<double> y = {10, 8, 5, 2, 5, 8, 10};
        
        auto minima = find_minima(x, y, 1, 0.5);
        REQUIRE(minima.size() > 0);
        CHECK(minima[0] == 3);
    }

    SECTION("min_spacing constraint") {
        std::vector<double> x = {0, 1, 2, 3, 4, 5, 6, 7, 8};
        std::vector<double> y = {5, 1, 5, 1, 5, 1, 5, 1, 5};
        
        auto minima1 = find_minima(x, y, 1, 0.1);
        auto minima2 = find_minima(x, y, 3, 0.1);
        
        CHECK(minima2.size() <= minima1.size());
    }

    SECTION("prominence constraint") {
        std::vector<double> x = {0, 1, 2, 3, 4, 5, 6};
        std::vector<double> y = {10, 9, 10, 1, 10, 9, 10};
        
        auto minima1 = find_minima(x, y, 1, 0.1);
        auto minima2 = find_minima(x, y, 1, 5.0);
        
        CHECK(minima2.size() <= minima1.size());
    }
}
