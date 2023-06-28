#include <catch2/catch_test_macros.hpp>

#include <utility/Limit.h>

TEST_CASE("Limit::Limit") {
    SECTION("default") {
        Limit limit;
        CHECK(limit.min == 0);
        CHECK(limit.max == 0);
    }

    SECTION("double, double") {
        Limit limit(1, 10);
        CHECK(limit.min == 1);
        CHECK(limit.max == 10);
    }
}

TEST_CASE("Limit::span") {
    Limit limit(1, 10);
    CHECK(limit.span() == 9);
}

TEST_CASE("Limit::center") {
    Limit limit(1, 10);
    CHECK(limit.center() == 5.5);
}

TEST_CASE("Limit::merge") {
    Limit limit1(1, 10);
    Limit limit2(2, 9);
    limit1.merge(limit2);
    CHECK(limit1.min == 1);
    CHECK(limit1.max == 10);

    Limit limit3(-1, 20);
    limit1.merge(limit3);
    CHECK(limit1.min == -1);
    CHECK(limit1.max == 20);
}

TEST_CASE("Limit::expand") {
    Limit limit(0, 10);
    limit.expand(0.5);
    CHECK(limit.min == -5);
    CHECK(limit.max == 15);
}

TEST_CASE("Limit::operator-=") {
    Limit limit(1, 10);
    limit -= 0.5;
    CHECK(limit.min == 0.5);
    CHECK(limit.max == 9.5);
}

TEST_CASE("Limit::operator-") {
    Limit limit1(1, 10);
    Limit limit2 = limit1 - 0.5;
    CHECK(limit2.min == 0.5);
    CHECK(limit2.max == 9.5);
}

TEST_CASE("Limit::operator+=") {
    Limit limit(1, 10);
    limit += 0.5;
    CHECK(limit.min == 1.5);
    CHECK(limit.max == 10.5);
}

TEST_CASE("Limit::operator+") {
    Limit limit1(1, 10);
    Limit limit2 = limit1 + 0.5;
    CHECK(limit2.min == 1.5);
    CHECK(limit2.max == 10.5);
}

TEST_CASE("Limit::operator==") {
    Limit limit1(1, 10);
    Limit limit2(1, 10);
    CHECK(limit1 == limit2);

    Limit limit3(1, 9);
    CHECK_FALSE(limit1 == limit3);
}

TEST_CASE("Limit::operator!=") {
    Limit limit1(1, 10);
    Limit limit2(1, 10);
    CHECK_FALSE(limit1 != limit2);

    Limit limit3(1, 9);
    CHECK(limit1 != limit3);
}

TEST_CASE("Limit::empty") {
    Limit limit1(1, 10);
    CHECK_FALSE(limit1.empty());

    Limit limit2(0, 0);
    CHECK(limit2.empty());
}