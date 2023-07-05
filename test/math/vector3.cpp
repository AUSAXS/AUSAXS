#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>

#include <vector>
#include <string>
#include <iostream>
#include <random>

#include <math/Vector3.h>

static Vector<double> GenRandVector(int m) {
    Vector<double> v(m);
    for (int i = 0; i < m; i++)
        v[i] = rand() % 100;
    return v;
}

TEST_CASE("basic operations") {
    Vector3<double> x = {1, 2, 3};
    Vector3<double> y = {4, 5, 6};
    Vector3<double> z = {7, 8, 9};

    // access
    REQUIRE(x.x() == 1);
    REQUIRE(x.y() == 2);
    REQUIRE(x.z() == 3);

    REQUIRE(y.x() == 4);
    REQUIRE(y.y() == 5);
    REQUIRE(y.z() == 6);

    // addition
    REQUIRE(x+y == Vector3{5, 7, 9});
    REQUIRE(x+z == Vector3{8, 10, 12});
    REQUIRE(y+z == Vector3{11, 13, 15});

    // subtraction
    REQUIRE(x-y == Vector3{-3, -3, -3});
    REQUIRE(x-z == Vector3{-6, -6, -6});
    REQUIRE(y-z == Vector3{-3, -3, -3});

    REQUIRE(y-x == Vector3{3, 3, 3});
    REQUIRE(z-x == Vector3{6, 6, 6});
    REQUIRE(z-y == Vector3{3, 3, 3});

    // division
    REQUIRE(x/2 == Vector3<double>{0.5, 1, 1.5});
    REQUIRE(y/2 == Vector3<double>{2, 2.5, 3});
    REQUIRE(z/2 == Vector3<double>{3.5, 4, 4.5});

    REQUIRE(2/x == Vector3<double>{2, 1, 2.0/3});
    REQUIRE(2/y == Vector3<double>{0.5, 2.0/5, 1.0/3});
    REQUIRE(2/z == Vector3<double>{2.0/7, 0.5, 2.0/9});

    // multiplication
    REQUIRE(x*2 == Vector3<double>{2, 4, 6});
    REQUIRE(y*2 == Vector3<double>{8, 10, 12});
    REQUIRE(z*2 == Vector3<double>{14, 16, 18});

    // negation
    REQUIRE(-x == Vector3{-1, -2, -3});
    REQUIRE(-y == Vector3{-4, -5, -6});
    REQUIRE(-z == Vector3{-7, -8, -9});

    // dot product
    REQUIRE(x.dot(y) == 4+10+18);
    REQUIRE(x.dot(z) == 7+16+27);
    REQUIRE(y.dot(z) == 28+40+54);

    // norm
    REQUIRE(x.norm() == sqrt(1+4+9));
    REQUIRE(y.norm() == sqrt(16+25+36));
    REQUIRE(z.norm() == sqrt(49+64+81));

    // normalize
    x.normalize(); y.normalize(); z.normalize();
    REQUIRE(x == Vector3<double>{0.2672612419124244, 0.5345224838248488, 0.8017837257372732});
    REQUIRE(y == Vector3<double>{0.4558423058385518, 0.5698028822981898, 0.6837634587578276});
    REQUIRE(z == Vector3<double>{0.5025707110324167, 0.5743665268941905, 0.6461623427559643});
}

TEST_CASE("constructors") {
    SECTION("vector3") {
        Vector3<double> x = {1, 2, 3};

        auto y = Vector3(x);
        REQUIRE(y == x);
        
        auto z = Vector3(std::move(x));
        REQUIRE(z == x);    

        Vector3<double> a;
        REQUIRE(a == Vector3{0, 0, 0});        
    }

    SECTION("vector") {
        Vector<double> x = {1, 2, 3};
        auto y = Vector3(x);
        REQUIRE(y == Vector3{1, 2, 3});

        auto z = Vector3(std::move(x));
        REQUIRE(z == Vector3{1, 2, 3});

        Vector<double> xx = {1, 2, 3, 4, 5};
        REQUIRE_THROWS(Vector3(xx));
    }

    SECTION("matrix") {
        Matrix<double> x = {
            {1, 2, 3}
        };
        auto y = Vector3(x);
        REQUIRE(y == Vector3{1, 2, 3});

        auto z = Vector3(std::move(x));
        REQUIRE(z == Vector3{1, 2, 3});

        Matrix<double> xx = {
            {1, 2, 3},
            {4, 5, 6}
        };
        REQUIRE_THROWS(Vector3(xx));
    }
}

TEST_CASE("assignment") {
    Vector3<double> x = {1, 2, 3};
    Vector3<double> y = {4, 5, 6};
    Vector3<double> z = {7, 8, 9};

    x += z;
    y += z;
    REQUIRE(x == Vector3<double>{8, 10, 12});
    REQUIRE(y == Vector3<double>{11, 13, 15});

    x -= z;
    y -= z;
    REQUIRE(x == Vector3<double>{1, 2, 3});
    REQUIRE(y == Vector3<double>{4, 5, 6});

    x = z;
    y = z;
    z = {0, 0, 0};
    REQUIRE(x == Vector3<double>{7, 8, 9});
    REQUIRE(y == Vector3<double>{7, 8, 9});
}

TEST_CASE("conversion") {
    Vector3<double> x = {1, 2, 3};

    Vector<double> v = x;
    REQUIRE(v == Vector<double>{1, 2, 3});
    REQUIRE(v.size() == 3);

    Matrix<double> M = x;
    REQUIRE(M == Matrix<double>{
        {1},
        {2}, 
        {3}
    });
    REQUIRE(M.N == 3);
    REQUIRE(M.M == 1);
}

TEST_CASE("matrix_product") {
    Vector3<double> x = {1, 2, 3};
    Vector3<double> y = {4, 5, 6};
    Vector3<double> z = {7, 8, 9};

    Matrix<double> m = {
        {1, 2, 3},
        {4, 5, 6},
        {7, 8, 9}
    };
    REQUIRE(m*x == Vector3<double>{14, 32, 50});
    REQUIRE(m*y == Vector3<double>{32, 77, 122});
    REQUIRE(m*z == Vector3<double>{50, 122, 194});
}

TEST_CASE("cross_product") {
    Vector3<double> x = {1, 2, 3};
    Vector3<double> y = {4, 5, 6};
    Vector3<double> z = {7, 8, 9};

    Vector3<double> v = x.cross(y);
    Vector3<double> w = x.cross(z);
    REQUIRE(v == Vector3<double>{-3, 6, -3});
    REQUIRE(w == Vector3<double>{-6, 12, -6});
}

TEST_CASE("distance") {
    Vector3<double> x = {1, 2, 3};
    Vector3<double> y = {4, 5, 6};
    Vector3<double> z = {7, 8, 9};

    REQUIRE(x.distance2(y) == 27);
    REQUIRE(x.distance2(z) == 108);
    REQUIRE(y.distance2(z) == 27);

    REQUIRE(x.distance(y) == sqrt(27));
    REQUIRE(x.distance(z) == sqrt(108));
    REQUIRE(y.distance(z) == sqrt(27));
}

TEST_CASE("rotation") {
    Vector3<double> x = {1, 2, 3};
    Vector3<double> y = {4, 5, 6};
    Vector3<double> z = {7, 8, 9};

    x = {1, 0, 0};
    y = {0, 1, 0};
    z = {0, 0, 1};

    Vector3<double> axis = {0, 1, 0};
    x.rotate(axis, M_PI_2);
    y.rotate(axis, M_PI_2);
    z.rotate(axis, M_PI_2);
    REQUIRE(x == Vector3<double>{0, 0, -1}); 
    REQUIRE(y == Vector3<double>{0, 1, 0}); 
    REQUIRE(z == Vector3<double>{1, 0, 0}); 

    axis = {1, 1, 1};
    x.rotate(axis, M_PI/4);
    y.rotate(axis, M_PI/4);
    z.rotate(axis, M_PI/4);
    REQUIRE(x == Vector3<double>{-0.5058793634, 0.3106172175, -0.8047378541}); 
    REQUIRE(y == Vector3<double>{-0.3106172175, 0.8047378541, 0.5058793634}); 
    REQUIRE(z == Vector3<double>{0.8047378541, 0.5058793634, -0.3106172175}); 

    x = {0, 2, 1};
    y = {5, 1, 3};
    z = {3, 7, 2};
    axis = {0.5, 2, 1};
    x.rotate(axis, 1.8);
    y.rotate(axis, 1.8);
    z.rotate(axis, 1.8);
    REQUIRE(x == Vector3<double>{0.5843819499, 1.6706126346, 1.3665837559}); 
    REQUIRE(y == Vector3<double>{1.8656722055, 4.7666664324, -2.9661689675}); 
    REQUIRE(z == Vector3<double>{0.0886646879, 7.4409765368, 2.5737145825}); 
}

TEST_CASE("normalize") {
    Vector3<double> x = {1, 2, 3};
    Vector3<double> y = {4, 5, 6};
    Vector3<double> z = {7, 8, 9};

    x = {2, 0, 0};
    REQUIRE(x.normalize() == Vector3<double>{1, 0, 0});

    x = {1, 1, 0};
    REQUIRE(x.normalize() == Vector3<double>{1, 1, 0}*sqrt(2)/2);
}

TEST_CASE("generate_basis") {
    Vector3<double> x = {1, 2, 3};
    Vector3<double> y = {4, 5, 6};
    Vector3<double> z = {7, 8, 9};

    x = {2, 0, 0};
    std::tie(x, y, z) = x.generate_basis();
    REQUIRE(x == Vector3<double>{1, 0, 0});
    REQUIRE((y == Vector3<double>{0, 1, 0} || y == Vector3<double>{0, 0, 1} || y == Vector3<double>{0, -1, 0} || y == Vector3<double>{0, 0, -1}));
    REQUIRE((z == Vector3<double>{0, 1, 0} || z == Vector3<double>{0, 0, 1} || z == Vector3<double>{0, -1, 0} || z == Vector3<double>{0, 0, -1}));

    for (int i = 0; i < 10; i++) {
        x = GenRandVector(3);
        std::tie(x, y, z) = x.generate_basis();
        REQUIRE_THAT(x.norm(), Catch::Matchers::WithinAbs(1, 1e-6));
        REQUIRE_THAT(y.norm(), Catch::Matchers::WithinAbs(1, 1e-6));
        REQUIRE_THAT(z.norm(), Catch::Matchers::WithinAbs(1, 1e-6));

        REQUIRE_THAT(x.dot(y), Catch::Matchers::WithinAbs(0, 1e-6));
        REQUIRE_THAT(y.dot(z), Catch::Matchers::WithinAbs(0, 1e-6));
        REQUIRE_THAT(z.dot(x), Catch::Matchers::WithinAbs(0, 1e-6));
    }
}