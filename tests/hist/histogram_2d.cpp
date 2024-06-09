#include <catch2/catch_test_macros.hpp>

#include <hist/Histogram2D.h>

TEST_CASE("Histogram2D::Histogram2D") {
    SECTION("default") {
        hist::Histogram2D hist;
        CHECK(hist.x_axis == Axis());
        CHECK(hist.y_axis == Axis());
    }

    SECTION("uint, uint") {
        hist::Histogram2D hist(10, 10);
        CHECK(hist.x_axis.bins == 10);
        CHECK(hist.y_axis.bins == 10);
    }

    SECTION("Axis&, Axis&") {
        Axis x_axis(1, 10, 2);
        Axis y_axis(1, 10, 5);
        hist::Histogram2D hist(x_axis, y_axis);
        CHECK(hist.x_axis == x_axis);
        CHECK(hist.y_axis == y_axis);
    }
}