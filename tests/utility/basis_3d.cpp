#include <catch2/catch_test_macros.hpp>

#include <utility/Basis3D.h>

using namespace ausaxs;

TEST_CASE("Basis3D::Basis3D") {
    SECTION("Vector3<double>&, Vector3<double>&, Vector3<double>&") {
        Vector3<double> x = {1, 2, 3};
        Vector3<double> y = {4, 5, 6};
        Vector3<double> z = {7, 8, 9};
        Basis3D basis(x, y, z);
        CHECK(basis.x == x);
        CHECK(basis.y == y);
        CHECK(basis.z == z);
    }
}