#include <catch2/catch_test_macros.hpp>

#include <utility/Axis.h>
#include <utility/Limit.h>

TEST_CASE("Axis::Axis") {
    SECTION("default") {
        Axis axis;
        CHECK(axis.span() == 0);
        CHECK(axis.empty());
    }

    SECTION("Limit&, int") {
        Limit limit(1, 10);
        Axis axis(limit, 3);
        CHECK(axis.width() == 3);
        CHECK(axis.span() == 9);
        CHECK(axis.step() == 3);
        CHECK_FALSE(axis.empty());
    }

    SECTION("int, double, double") {
        Axis axis(1, 10, 3);
        CHECK(axis.width() == 3);
        CHECK(axis.span() == 9);
        CHECK(axis.step() == 3);
        CHECK_FALSE(axis.empty());
    }
}

TEST_CASE("Axis::operator=") {
    Axis axis1(1, 9, 3);
    Axis axis2;
    axis2 = axis1;
    CHECK(axis1 == axis2);
}

TEST_CASE("Axis::operator==") {
    Axis axis1(1, 9, 3);
    Axis axis2(1, 9, 3);
    CHECK(axis1 == axis2);

    Axis axis3(1, 9, 4);
    CHECK_FALSE(axis1 == axis3);
}

TEST_CASE("Axis::operator!=") {
    Axis axis1(1, 9, 3);
    Axis axis2(1, 9, 3);
    CHECK_FALSE(axis1 != axis2);

    Axis axis3(1, 9, 4);
    CHECK(axis1 != axis3);
}

TEST_CASE("Axis::width") {
    Axis axis(1, 10, 3);
    CHECK(axis.width() == 3);
}

TEST_CASE("Axis::span") {
    Axis axis(1, 9, 3);
    CHECK(axis.span() == 8);
}

TEST_CASE("Axis::step") {
    Axis axis(1, 10, 3);
    CHECK(axis.step() == 3);
}

TEST_CASE("Axis::resize") {}

TEST_CASE("Axis::as_vector") {}

TEST_CASE("Axis::empty") {}

TEST_CASE("Axis::limits") {}