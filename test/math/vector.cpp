#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>

#include <vector>
#include <string>
#include <iostream>
#include <random>

#include <math/Vector.h>

static Vector<double> GenRandVector(int m) {
    Vector<double> v(m);
    for (int i = 0; i < m; i++)
        v[i] = rand() % 100;
    return v;
}

TEST_CASE("basic operations") {
    Vector<double> x = {1, 2, 3, 4};
    Vector<double> y = {2, 3, 4, 5};
    Vector<double> z = {4, 3, 2, 1};

    // addition
    REQUIRE(x+y == Vector<double>{3, 5, 7, 9});
    REQUIRE(x+z == Vector<double>{5, 5, 5, 5});
    REQUIRE(y+z == Vector<double>{6, 6, 6, 6});

    // subtraction
    REQUIRE(x-y == Vector<double>{-1, -1, -1, -1});
    REQUIRE(x-z == Vector<double>{-3, -1, 1, 3});
    REQUIRE(y-z == Vector<double>{-2, 0, 2, 4});

    REQUIRE(y-x == Vector<double>{1, 1, 1, 1});
    REQUIRE(z-x == Vector<double>{3, 1, -1, -3});
    REQUIRE(z-y == Vector<double>{2, 0, -2, -4});

    // negation
    REQUIRE(-x == Vector<double>{-1, -2, -3, -4});
    REQUIRE(-y == Vector<double>{-2, -3, -4, -5});
    REQUIRE(-z == Vector<double>{-4, -3, -2, -1});

    // dot product
    REQUIRE(x.dot(y) == 2+6+12+20);
    REQUIRE(x.dot(z) == 4+6+6+4);
    REQUIRE(y.dot(z) == 8+9+8+5);

    // norm
    REQUIRE(x.norm() == sqrt(1+4+9+16));
    REQUIRE(y.norm() == sqrt(4+9+16+25));
    REQUIRE(z.norm() == sqrt(16+9+4+1));

    // vector multiplication
    REQUIRE(x*y == Vector<double>{2, 6, 12, 20});
    REQUIRE(x*z == Vector<double>{4, 6, 6, 4});
    REQUIRE(y*z == Vector<double>{8, 9, 8, 5});

    // scalar multiplication
    REQUIRE(x*2 == Vector<double>{2, 4, 6, 8});
    REQUIRE(y*3 == Vector<double>{6, 9, 12, 15});
    REQUIRE(z*5 == Vector<double>{20, 15, 10, 5});

    REQUIRE(2*x == Vector<double>{2, 4, 6, 8});
    REQUIRE(3*y == Vector<double>{6, 9, 12, 15});
    REQUIRE(5*z == Vector<double>{20, 15, 10, 5});

    // scalar division
    REQUIRE(x/2 == Vector<double>{1./2, 2./2, 3./2, 4./2});
    REQUIRE(y/4 == Vector<double>{2./4, 3./4, 4./4, 5./4});
    REQUIRE(z/8 == Vector<double>{4./8, 3./8, 2./8, 1./8});
}

TEST_CASE("assignment") {
    Vector<double> x = {1, 2, 3, 4};
    Vector<double> y = {2, 3, 4, 5};
    Vector<double> z = {4, 3, 2, 1};

    x += z;
    y += z;
    REQUIRE(x == Vector<double>{5, 5, 5, 5});
    REQUIRE(y == Vector<double>{6, 6, 6, 6});

    x -= z;
    y -= z;
    REQUIRE(x == Vector<double>{1, 2, 3, 4});
    REQUIRE(y == Vector<double>{2, 3, 4, 5});

    x = z;
    y = z;
    z = {0, 0, 0};
    REQUIRE(x == Vector<double>{4, 3, 2, 1});
    REQUIRE(y == Vector<double>{4, 3, 2, 1});

    x = {9, 8, 7, 6, 5};
    y = {1, 2, 3, 4, 5};
    REQUIRE(x == Vector<double>{9, 8, 7, 6, 5});
    REQUIRE(y == Vector<double>{1, 2, 3, 4, 5});
}

TEST_CASE("distance") {
    Vector<double> x = {1, 2, 3, 4};
    Vector<double> y = {2, 3, 4, 5};
    Vector<double> z = {4, 3, 2, 1};

    REQUIRE(x.distance2(y) == 1+1+1+1);
    REQUIRE(x.distance2(z) == 9+1+1+9);
    REQUIRE(y.distance2(z) == 4+0+4+16);

    REQUIRE(x.distance(y) == sqrt(1+1+1+1));
    REQUIRE(x.distance(z) == sqrt(9+1+1+9));
    REQUIRE(y.distance(z) == sqrt(4+0+4+16));
}

TEST_CASE("iterator") {
    Vector<double> x = {1, 2, 3, 4};
    Vector<double> y = {2, 3, 4, 5};
    Vector<double> z = {4, 3, 2, 1};

    for (const auto& e : x) {
        REQUIRE((e == 1 || e == 2 || e == 3 || e == 4));
    }
    for (const auto& e : y) {
        REQUIRE((e == 2 || e == 3 || e == 4 || e == 5));
    }
}