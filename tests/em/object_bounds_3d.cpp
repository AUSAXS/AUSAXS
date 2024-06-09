#include <catch2/catch_test_macros.hpp>

#include <em/ObjectBounds3D.h>
#include <em/ObjectBounds2D.h>

TEST_CASE("ObjectBounds3D::ObjectBounds3D") {
    SECTION("unsigned int, unsigned int, unsigned int") {
        em::ObjectBounds3D bounds(1, 2, 3);
        CHECK(bounds.size_x() == 1);
        CHECK(bounds.size_y() == 2);
        CHECK(bounds.size_z() == 3);
    }
}

TEST_CASE("ObjectBounds3D::operator[]") {
    em::ObjectBounds3D bounds(1, 2, 3);
    CHECK(bounds[0].size_x() == 1);
    CHECK(bounds[0].size_y() == 2);
}

TEST_CASE("ObjectBounds3D::bounded_volume") {
    em::ObjectBounds3D bounds(1, 2, 3);
    CHECK(bounds.bounded_volume() == 9);

    bounds[0].set_bounds(0, 0, 1);
    CHECK(bounds.bounded_volume() == 8);

    bounds[0].set_bounds(0, 0, 0);
    CHECK(bounds.bounded_volume() == 7);
}

TEST_CASE("ObjectBounds3D::total_volume") {
    em::ObjectBounds3D bounds(1, 2, 3);
    CHECK(bounds.total_volume() == 6);

    bounds[0].set_bounds(0, 0, 1);
    CHECK(bounds.total_volume() == 6);

    bounds[0].set_bounds(0, 0, 0);
    CHECK(bounds.total_volume() == 6);
}