#include <catch2/catch_test_macros.hpp>

#include <utility/Axis3D.h>
#include <utility/Limit3D.h>
#include <math/Vector3.h>

TEST_CASE("Axis3D::Axis3D") {
    SECTION("default") {
        Axis3D axis;
        CHECK(axis.x == Axis());
        CHECK(axis.y == Axis());
        CHECK(axis.z == Axis());
    }

    SECTION("Axis3D&") {
        Axis3D axis1;
        Axis3D axis2(axis1);
        CHECK(axis1.x == axis2.x);
    }

    SECTION("Axis&, Axis&, Axis&") {
        Axis x(0, 1, 1);
        Axis y(2, 5, 2);
        Axis z(4, 9, 3);
        Axis3D axis(x, y, z);
        CHECK(axis.x == x);
        CHECK(axis.y == y);
        CHECK(axis.z == z);
    }

    SECTION("Limit3D&, double") {
        Limit3D limits(0, 1, 2, 3, 4, 5);
        double width = 1;
        Axis3D axis(limits, width);
        CHECK(axis.x == Axis(0, 1, 1));
        CHECK(axis.y == Axis(2, 3, 1));
        CHECK(axis.z == Axis(4, 5, 1));
    }

    SECTION("double, double, double, double, double, double, double") {
        double xmin = 0;
        double xmax = 1;
        double ymin = 2;
        double ymax = 3;
        double zmin = 4;
        double zmax = 5;
        double width = 1;
        Axis3D axis(xmin, xmax, ymin, ymax, zmin, zmax, width);
        CHECK(axis.x == Axis(xmin, xmax, 1));
        CHECK(axis.y == Axis(ymin, ymax, 1));
        CHECK(axis.z == Axis(zmin, zmax, 1));
    }

    SECTION("Vector3<double>&, Vector3<double>&, double") {
        Vector3<double> min = {0, 2, 4};
        Vector3<double> max = {1, 3, 5};
        double width = 1;
        Axis3D axis(min, max, width);
        CHECK(axis.x == Axis(0, 1, 1));
        CHECK(axis.y == Axis(2, 3, 1));
        CHECK(axis.z == Axis(4, 5, 1));
    }
}

TEST_CASE("Axis3D::operator=") {
    Axis3D axis1(0, 1, 2, 3, 4, 5, 1);
    Axis3D axis2;
    axis2 = axis1;
    CHECK(axis1.x == axis2.x);
}

TEST_CASE("Axis3D::operator==") {
    Axis3D axis1(0, 1, 2, 3, 4, 5, 1);
    Axis3D axis2(0, 1, 2, 3, 4, 5, 1);
    CHECK(axis1 == axis2);

    Axis3D axis3(0, 1, 2, 3, 4, 6, 1);
    CHECK_FALSE(axis1 == axis3);
}

TEST_CASE("Axis3D::operator!=") {
    Axis3D axis1(0, 1, 2, 3, 4, 5, 1);
    Axis3D axis2(0, 1, 2, 3, 4, 5, 1);
    CHECK_FALSE(axis1 != axis2);

    Axis3D axis3(0, 1, 2, 3, 4, 6, 1);
    CHECK(axis1 != axis3);
}

TEST_CASE("Axis3D::empty") {
    Axis3D axis;
    CHECK(axis.empty());

    Axis3D axis2(0, 1, 2, 3, 4, 5, 1);
    CHECK_FALSE(axis2.empty());
}

TEST_CASE("Axis3D::rebin") {
    Axis3D axis(0, 2, 2, 6, 4, 10, 1);
    axis.rebin(2);
    CHECK(axis.x == Axis(0, 2, 1));
    CHECK(axis.y == Axis(2, 6, 2));
    CHECK(axis.z == Axis(4, 10, 3));
}

TEST_CASE("Axis3D::width") {
    Axis3D axis(0, 1, 2, 3, 4, 5, 1);
    CHECK(axis.width() == 1);

    Axis3D axis2(0, 2, 2, 6, 4, 10, 2);
    CHECK(axis2.width() == 2);

    axis2.x = Axis(0, 5, 3);
    CHECK_THROWS(axis2.width()); // inconsistent widths
}