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

    SECTION("single point minimum") {
        std::vector<double> x = {0, 1, 2, 3, 4, 5, 6};
        std::vector<double> y = {10, 8, 5, 2, 5, 8, 10};
        
        auto minima = find_minima(x, y, 1, 0.5);
        REQUIRE(minima.size() > 0);
        CHECK(minima[0] == 3);
    }
}
