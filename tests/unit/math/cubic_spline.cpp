#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <math/CubicSpline.h>

using namespace ausaxs::math;

TEST_CASE("CubicSpline::CubicSpline") {
    SECTION("linear data") {
        std::vector<double> x = {0, 1, 2, 3, 4};
        std::vector<double> y = {0, 1, 2, 3, 4};
        
        CubicSpline spline(x, y);
        
        CHECK_THAT(spline.spline(0), Catch::Matchers::WithinAbs(0, 1e-6));
        CHECK_THAT(spline.spline(1), Catch::Matchers::WithinAbs(1, 1e-6));
        CHECK_THAT(spline.spline(2), Catch::Matchers::WithinAbs(2, 1e-6));
        CHECK_THAT(spline.spline(3), Catch::Matchers::WithinAbs(3, 1e-6));
        CHECK_THAT(spline.spline(4), Catch::Matchers::WithinAbs(4, 1e-6));
    }

    SECTION("quadratic data") {
        std::vector<double> x = {0, 1, 2, 3};
        std::vector<double> y = {0, 1, 4, 9};
        
        CubicSpline spline(x, y);
        
        CHECK_THAT(spline.spline(0), Catch::Matchers::WithinAbs(0, 1e-3));
        CHECK_THAT(spline.spline(1), Catch::Matchers::WithinAbs(1, 1e-3));
        CHECK_THAT(spline.spline(2), Catch::Matchers::WithinAbs(4, 1e-3));
        CHECK_THAT(spline.spline(3), Catch::Matchers::WithinAbs(9, 1e-3));
    }

    SECTION("interpolation") {
        std::vector<double> x = {0, 1, 2, 3, 4};
        std::vector<double> y = {0, 1, 4, 9, 16};
        
        CubicSpline spline(x, y);
        
        double mid = spline.spline(0.5);
        CHECK(mid > 0);
        CHECK(mid < 1);
    }
}

TEST_CASE("CubicSpline::spline") {
    SECTION("on grid points") {
        std::vector<double> x = {1, 2, 3, 4, 5};
        std::vector<double> y = {2, 4, 6, 8, 10};
        
        CubicSpline spline(x, y);
        
        for (unsigned int i = 0; i < x.size(); ++i) {
            CHECK_THAT(spline.spline(x[i]), Catch::Matchers::WithinAbs(y[i], 1e-6));
        }
    }

    SECTION("between grid points") {
        std::vector<double> x = {0, 1, 2, 3};
        std::vector<double> y = {0, 1, 0, 1};
        
        CubicSpline spline(x, y);
        
        double val = spline.spline(1.5);
        CHECK(val >= 0);
        CHECK(val <= 1);
    }
}
