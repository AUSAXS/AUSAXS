#include <catch2/catch_test_macros.hpp>

#include <utility/Limit3D.h>

TEST_CASE("Limit3D::Limit3D") {
    SECTION("default") {
        Limit3D limit;
        CHECK(limit.x.min == 0);
        CHECK(limit.x.max == 0);
        CHECK(limit.y.min == 0);
        CHECK(limit.y.max == 0);
        CHECK(limit.z.min == 0);
        CHECK(limit.z.max == 0);
    }

    SECTION("Limit&, Limit&, Limit&") {
        Limit x(0, 1);
        Limit y(2, 3);
        Limit z(4, 5);
        Limit3D limit(x, y, z);
        CHECK(limit.x == x);
        CHECK(limit.y == y);
        CHECK(limit.z == z);
    }

    SECTION("double, double, double, double, double, double") {
        double xmin = 0;
        double xmax = 1;
        double ymin = 2;
        double ymax = 3;
        double zmin = 4;
        double zmax = 5;
        Limit3D limit(xmin, xmax, ymin, ymax, zmin, zmax);
        CHECK(limit.x == Limit(xmin, xmax));
        CHECK(limit.y == Limit(ymin, ymax));
        CHECK(limit.z == Limit(zmin, zmax));
    }
}        

TEST_CASE("Limit3D::empty") {
    Limit3D limit;
    CHECK(limit.empty() == true);

    Limit3D limit2(1, 2, 3, 4, 5, 6);
    CHECK(limit2.empty() == false);
}