#include <catch2/catch_test_macros.hpp>

#include <em/ObjectBounds2D.h>
#include <utility/Limit.h>

TEST_CASE("ObjectBounds2D::ObjectBounds2D") {
    SECTION("unsigned int, unsigned int") {
        em::ObjectBounds2D bounds = em::ObjectBounds2D(1, 2);
        REQUIRE(bounds.size_x() == 1);
        REQUIRE(bounds.size_y() == 2);
    }
}

TEST_CASE("ObjectBounds2D::operator[]") {
    em::ObjectBounds2D bounds = em::ObjectBounds2D(2, 2);
    bounds.set_bounds(0, Limit(0, 1));
    bounds.set_bounds(1, Limit(1, 1));
    CHECK(bounds[0] == Limit(0, 1));
    CHECK(bounds[1] == Limit(1, 1));
}

TEST_CASE("ObjectBounds2D::set_min") {
    em::ObjectBounds2D bounds = em::ObjectBounds2D(2, 2);
    bounds.set_min(0, 1);
    CHECK(bounds[0] == Limit(1, 2));

    CHECK_THROWS(bounds.set_min(2, 0));
    CHECK_THROWS(bounds.set_min(0, 3));
}

TEST_CASE("ObjectBounds2D::set_max") {
    em::ObjectBounds2D bounds = em::ObjectBounds2D(2, 2);
    bounds.set_max(0, 1);
    CHECK(bounds[0] == Limit(0, 1));

    CHECK_THROWS(bounds.set_max(2, 0));
    CHECK_THROWS(bounds.set_max(0, 3));
}

TEST_CASE("ObjectBounds2D::set_bounds") {
    SECTION("unsigned int, Limit") {
        em::ObjectBounds2D bounds = em::ObjectBounds2D(2, 2);
        bounds.set_bounds(0, Limit(0, 1));
        bounds.set_bounds(1, Limit(1, 1));
        CHECK(bounds[0] == Limit(0, 1));
        CHECK(bounds[1] == Limit(1, 1));

        bounds.set_bounds(1, Limit(0, 2));
        CHECK(bounds[1] == Limit(0, 2));

        CHECK_THROWS(bounds.set_bounds(2, Limit(0, 1)));
        CHECK_THROWS(bounds.set_bounds(0, Limit(1, 3)));
    }

    SECTION("unsigned int, unsigned int, unsigned int") {
        em::ObjectBounds2D bounds = em::ObjectBounds2D(2, 2);
        bounds.set_bounds(0, 0, 1);
        bounds.set_bounds(1, 1, 1);
        CHECK(bounds[0] == Limit(0, 1));
        CHECK(bounds[1] == Limit(1, 1));

        bounds.set_bounds(1, 0, 2);
        CHECK(bounds[1] == Limit(0, 2));

        CHECK_THROWS(bounds.set_bounds(2, 0, 1));
        CHECK_THROWS(bounds.set_bounds(0, 1, 3));
    }
}

TEST_CASE("ObjectBounds2D::empty") {
    em::ObjectBounds2D bounds = em::ObjectBounds2D(2, 2);
    CHECK(!bounds.empty());    
    bounds.set_bounds(0, Limit(0, 0));
    bounds.set_bounds(1, Limit(0, 0));
    CHECK(bounds.empty());
}

TEST_CASE("ObjectBounds2D::bounded_area") {
    em::ObjectBounds2D bounds = em::ObjectBounds2D(2, 5);
    CHECK(bounds.bounded_area() == 12);

    bounds.set_bounds(0, Limit(0, 1));
    bounds.set_bounds(1, Limit(1, 1));
    CHECK(bounds.bounded_area() == 3);

    bounds.set_bounds(0, Limit(0, 5));
    bounds.set_bounds(1, Limit(3, 5));
    CHECK(bounds.bounded_area() == 9);    
}

TEST_CASE("ObjectBounds2D::total_area") {
    em::ObjectBounds2D bounds = em::ObjectBounds2D(2, 5);
    CHECK(bounds.total_area() == 10);

    bounds.set_bounds(0, Limit(0, 1));
    bounds.set_bounds(1, Limit(1, 1));
    CHECK(bounds.total_area() == 10);

    bounds.set_bounds(0, Limit(0, 5));
    bounds.set_bounds(1, Limit(3, 5));
    CHECK(bounds.total_area() == 10);
}

TEST_CASE("ObjectBounds2D::operator==") {
    em::ObjectBounds2D bounds1 = em::ObjectBounds2D(2, 5);
    em::ObjectBounds2D bounds2 = em::ObjectBounds2D(2, 5);
    CHECK(bounds1 == bounds2);

    bounds1.set_bounds(0, Limit(0, 1));
    bounds1.set_bounds(1, Limit(1, 1));
    CHECK(bounds1 != bounds2);

    bounds2.set_bounds(0, Limit(0, 1));
    bounds2.set_bounds(1, Limit(1, 1));
    CHECK(bounds1 == bounds2);

    bounds1.set_bounds(0, Limit(0, 5));
    bounds1.set_bounds(1, Limit(3, 5));
    CHECK(bounds1 != bounds2);

    bounds2.set_bounds(0, Limit(0, 5));
    bounds2.set_bounds(1, Limit(3, 5));
    CHECK(bounds1 == bounds2);
}