#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <dataset/PointSet.h>

using namespace ausaxs;

TEST_CASE("Point1D::Point1D") {
    SECTION("default constructor") {
        Point1D p;
        CHECK(p.x == 0);
        CHECK(p.xerr == 0);
    }

    SECTION("double") {
        Point1D p(5.0);
        CHECK(p.x == 5.0);
        CHECK(p.xerr == 0);
    }

    SECTION("double, double") {
        Point1D p(5.0, 1.0);
        CHECK(p.x == 5.0);
        CHECK(p.xerr == 1.0);
    }
}

TEST_CASE("Point1D::dim") {
    CHECK(Point1D::dim() == 1);
}

TEST_CASE("Point1D::operator==") {
    Point1D p1(5.0, 1.0);
    Point1D p2(5.0, 1.0);
    Point1D p3(5.0, 2.0);
    Point1D p4(6.0, 1.0);

    CHECK(p1 == p2);
    CHECK_FALSE(p1 == p3);
    CHECK_FALSE(p1 == p4);
}

TEST_CASE("Point1D::operator!=") {
    Point1D p1(5.0, 1.0);
    Point1D p2(5.0, 1.0);
    Point1D p3(5.0, 2.0);

    CHECK_FALSE(p1 != p2);
    CHECK(p1 != p3);
}

TEST_CASE("Point2D::Point2D") {
    SECTION("default constructor") {
        Point2D p;
        CHECK(p.x == 0);
        CHECK(p.y == 0);
        CHECK(p.xerr == 0);
        CHECK(p.yerr == 0);
    }

    SECTION("double, double") {
        Point2D p(5.0, 10.0);
        CHECK(p.x == 5.0);
        CHECK(p.y == 10.0);
        CHECK(p.xerr == 0);
        CHECK(p.yerr == 0);
    }

    SECTION("double, double, double") {
        Point2D p(5.0, 10.0, 2.0);
        CHECK(p.x == 5.0);
        CHECK(p.y == 10.0);
        CHECK(p.xerr == 0);
        CHECK(p.yerr == 2.0);
    }

    SECTION("double, double, double, double") {
        Point2D p(5.0, 10.0, 1.0, 2.0);
        CHECK(p.x == 5.0);
        CHECK(p.y == 10.0);
        CHECK(p.xerr == 1.0);
        CHECK(p.yerr == 2.0);
    }
}

TEST_CASE("Point2D::dim") {
    CHECK(Point2D::dim() == 2);
}

TEST_CASE("Point2D::operator==") {
    Point2D p1(5.0, 10.0, 1.0, 2.0);
    Point2D p2(5.0, 10.0, 1.0, 2.0);
    Point2D p3(5.0, 10.0, 1.0, 3.0);
    Point2D p4(5.0, 11.0, 1.0, 2.0);
    Point2D p5(6.0, 10.0, 1.0, 2.0);
    Point2D p6(5.0, 10.0, 2.0, 2.0);

    CHECK(p1 == p2);
    CHECK_FALSE(p1 == p3);
    CHECK_FALSE(p1 == p4);
    CHECK_FALSE(p1 == p5);
    CHECK_FALSE(p1 == p6);
}

TEST_CASE("Point2D::operator!=") {
    Point2D p1(5.0, 10.0, 1.0, 2.0);
    Point2D p2(5.0, 10.0, 1.0, 2.0);
    Point2D p3(5.0, 10.0, 1.0, 3.0);

    CHECK_FALSE(p1 != p2);
    CHECK(p1 != p3);
}

TEST_CASE("Point3D::Point3D") {
    SECTION("default constructor") {
        Point3D p;
        CHECK(p.x == 0);
        CHECK(p.y == 0);
        CHECK(p.z == 0);
        CHECK(p.xerr == 0);
        CHECK(p.yerr == 0);
        CHECK(p.zerr == 0);
    }

    SECTION("double, double, double") {
        Point3D p(5.0, 10.0, 15.0);
        CHECK(p.x == 5.0);
        CHECK(p.y == 10.0);
        CHECK(p.z == 15.0);
        CHECK(p.xerr == 0);
        CHECK(p.yerr == 0);
        CHECK(p.zerr == 0);
    }
}

TEST_CASE("Point3D::dim") {
    CHECK(Point3D::dim() == 3);
}

TEST_CASE("Point3D::operator==") {
    Point3D p1(5.0, 10.0, 15.0);
    Point3D p2(5.0, 10.0, 15.0);
    Point3D p3(5.0, 10.0, 16.0);
    Point3D p4(5.0, 11.0, 15.0);
    Point3D p5(6.0, 10.0, 15.0);

    CHECK(p1 == p2);
    CHECK_FALSE(p1 == p3);
    CHECK_FALSE(p1 == p4);
    CHECK_FALSE(p1 == p5);
}

TEST_CASE("Point3D::operator!=") {
    Point3D p1(5.0, 10.0, 15.0);
    Point3D p2(5.0, 10.0, 15.0);
    Point3D p3(5.0, 10.0, 16.0);

    CHECK_FALSE(p1 != p2);
    CHECK(p1 != p3);
}
