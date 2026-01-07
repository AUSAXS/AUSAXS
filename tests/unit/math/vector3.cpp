#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <math/Vector3.h>
#include <math/Vector.h>
#include <math/Matrix.h>

#include <numbers>

using namespace ausaxs;

TEST_CASE("Vector3::Vector3") {
    SECTION("xyz constructor") {
        Vector3<double> x(1, 2, 3);
        REQUIRE(x.x() == 1);
        REQUIRE(x.y() == 2);
        REQUIRE(x.z() == 3);
    }

    SECTION("initializer_list") {
        Vector3<double> x = {1, 2, 3};
        REQUIRE(x.x() == 1);
        REQUIRE(x.y() == 2);
        REQUIRE(x.z() == 3);
    }

    SECTION("copy constructor") {
        Vector3<double> x = {1, 2, 3};
        Vector3<double> y(x);
        REQUIRE(y == x);
    }

    SECTION("move constructor") {
        Vector3<double> x = {1, 2, 3};
        Vector3<double> y(std::move(x));
        REQUIRE(y.x() == 1);
        REQUIRE(y.y() == 2);
        REQUIRE(y.z() == 3);
    }

    SECTION("from Vector") {
        Vector<double> v = {1, 2, 3};
        Vector3<double> x(v);
        REQUIRE(x.x() == 1);
        REQUIRE(x.y() == 2);
        REQUIRE(x.z() == 3);
    }

    SECTION("from Matrix") {
        Matrix<double> M = {{1}, {2}, {3}};
        Vector3<double> x(M);
        REQUIRE(x.x() == 1);
        REQUIRE(x.y() == 2);
        REQUIRE(x.z() == 3);
    }
}

TEST_CASE("Vector3::operator[]") {
    Vector3<double> x = {1, 2, 3};
    
    SECTION("const access") {
        const Vector3<double>& cx = x;
        REQUIRE(cx[0] == 1);
        REQUIRE(cx[1] == 2);
        REQUIRE(cx[2] == 3);
    }

    SECTION("mutable access") {
        x[0] = 10;
        x[1] = 20;
        x[2] = 30;
        REQUIRE(x[0] == 10);
        REQUIRE(x[1] == 20);
        REQUIRE(x[2] == 30);
    }
}

TEST_CASE("Vector3::x/y/z") {
    Vector3<double> v = {1, 2, 3};
    REQUIRE(v.x() == 1);
    REQUIRE(v.y() == 2);
    REQUIRE(v.z() == 3);
}

TEST_CASE("Vector3::operator=") {
    SECTION("copy assignment") {
        Vector3<double> x = {1, 2, 3};
        Vector3<double> y;
        y = x;
        REQUIRE(y == x);
    }

    SECTION("move assignment") {
        Vector3<double> x = {1, 2, 3};
        Vector3<double> y;
        y = std::move(x);
        REQUIRE(y.x() == 1);
        REQUIRE(y.y() == 2);
        REQUIRE(y.z() == 3);
    }

    SECTION("initializer_list assignment") {
        Vector3<double> x;
        x = {1, 2, 3};
        REQUIRE(x.x() == 1);
        REQUIRE(x.y() == 2);
        REQUIRE(x.z() == 3);
    }
}

TEST_CASE("Vector3::operator==") {
    Vector3<double> x = {1, 2, 3};
    Vector3<double> y = {1, 2, 3};
    Vector3<double> z = {1, 2, 4};

    REQUIRE(x == y);
    REQUIRE_FALSE(x == z);
}

TEST_CASE("Vector3::operator!=") {
    Vector3<double> x = {1, 2, 3};
    Vector3<double> y = {1, 2, 3};
    Vector3<double> z = {1, 2, 4};

    REQUIRE_FALSE(x != y);
    REQUIRE(x != z);
}

TEST_CASE("Vector3::operator+") {
    Vector3<double> x = {1, 2, 3};
    Vector3<double> y = {4, 5, 6};
    
    Vector3<double> result = x + y;
    REQUIRE(result == Vector3<double>{5, 7, 9});
}

TEST_CASE("Vector3::operator-") {
    Vector3<double> x = {1, 2, 3};
    Vector3<double> y = {4, 5, 6};
    
    Vector3<double> result = x - y;
    REQUIRE(result == Vector3<double>{-3, -3, -3});
}

TEST_CASE("Vector3::operator- (unary)") {
    Vector3<double> x = {1, 2, 3};
    Vector3<double> result = -x;
    REQUIRE(result == Vector3<double>{-1, -2, -3});
}

TEST_CASE("Vector3::operator*") {
    Vector3<double> x = {1, 2, 3};
    
    SECTION("right multiplication") {
        Vector3<double> result = x * 2;
        REQUIRE(result == Vector3<double>{2, 4, 6});
    }

    SECTION("left multiplication") {
        Vector3<double> result = 2 * x;
        REQUIRE(result == Vector3<double>{2, 4, 6});
    }
}

TEST_CASE("Vector3::operator/") {
    Vector3<double> x = {2, 4, 6};
    
    SECTION("right division") {
        Vector3<double> result = x / 2;
        REQUIRE(result == Vector3<double>{1, 2, 3});
    }

    SECTION("left division") {
        Vector3<double> result = 6 / x;
        REQUIRE(result == Vector3<double>{3, 1.5, 1});
    }
}

TEST_CASE("Vector3::operator+=") {
    Vector3<double> x = {1, 2, 3};
    Vector3<double> y = {4, 5, 6};
    
    x += y;
    REQUIRE(x == Vector3<double>{5, 7, 9});
}

TEST_CASE("Vector3::operator-=") {
    Vector3<double> x = {1, 2, 3};
    Vector3<double> y = {4, 5, 6};
    
    x -= y;
    REQUIRE(x == Vector3<double>{-3, -3, -3});
}

TEST_CASE("Vector3::operator*=") {
    Vector3<double> x = {1, 2, 3};
    x *= 2;
    REQUIRE(x == Vector3<double>{2, 4, 6});
}

TEST_CASE("Vector3::operator/=") {
    Vector3<double> x = {2, 4, 6};
    x /= 2;
    REQUIRE(x == Vector3<double>{1, 2, 3});
}

TEST_CASE("Vector3::dot") {
    Vector3<double> x = {1, 2, 3};
    Vector3<double> y = {4, 5, 6};
    
    double result = x.dot(y);
    REQUIRE(result == 4+10+18);
}

TEST_CASE("Vector3::cross") {
    Vector3<double> x = {1, 2, 3};
    Vector3<double> y = {4, 5, 6};
    
    Vector3<double> result = x.cross(y);
    REQUIRE(result == Vector3<double>{-3, 6, -3});
}

TEST_CASE("Vector3::norm") {
    Vector3<double> x = {1, 2, 3};
    double result = x.norm();
    REQUIRE(result == sqrt(1+4+9));
}

TEST_CASE("Vector3::distance") {
    Vector3<double> x = {1, 2, 3};
    Vector3<double> y = {4, 5, 6};
    
    double result = x.distance(y);
    REQUIRE(result == sqrt(27));
}

TEST_CASE("Vector3::distance2") {
    Vector3<double> x = {1, 2, 3};
    Vector3<double> y = {4, 5, 6};
    
    double result = x.distance2(y);
    REQUIRE(result == 27);
}

TEST_CASE("Vector3::normalize") {
    SECTION("non-zero vector") {
        Vector3<double> x = {2, 0, 0};
        Vector3<double> result = x.normalize();
        REQUIRE(result == Vector3<double>{1, 0, 0});
        REQUIRE_THAT(result.norm(), Catch::Matchers::WithinAbs(1.0, 1e-10));
    }

    SECTION("arbitrary vector") {
        Vector3<double> x = {1, 1, 0};
        Vector3<double> result = x.normalize();
        REQUIRE_THAT(result.norm(), Catch::Matchers::WithinAbs(1.0, 1e-10));
        REQUIRE(result == Vector3<double>{1, 1, 0}*sqrt(2)/2);
    }
}

TEST_CASE("Vector3::rotate") {
    SECTION("90 degree rotation around y-axis") {
        Vector3<double> x = {1, 0, 0};
        Vector3<double> axis = {0, 1, 0};
        x.rotate(axis, std::numbers::pi/2);
        REQUIRE(x == Vector3<double>{0, 0, -1}); 
    }

    SECTION("identity rotation") {
        Vector3<double> x = {1, 2, 3};
        Vector3<double> axis = {0, 1, 0};
        x.rotate(axis, 0);
        REQUIRE(x == Vector3<double>{1, 2, 3});
    }
}

TEST_CASE("Vector3::conversion to Vector") {
    Vector3<double> x = {1, 2, 3};
    Vector<double> v = x;
    REQUIRE(v.size() == 3);
    REQUIRE(v[0] == 1);
    REQUIRE(v[1] == 2);
    REQUIRE(v[2] == 3);
}

TEST_CASE("Vector3::conversion to Matrix") {
    Vector3<double> x = {1, 2, 3};
    Matrix<double> M = x;
    REQUIRE(M.N == 3);
    REQUIRE(M.M == 1);
    REQUIRE(M[0][0] == 1);
    REQUIRE(M[1][0] == 2);
    REQUIRE(M[2][0] == 3);
}

TEST_CASE("Vector3::structured bindings") {
    Vector3<double> x = {1, 2, 3};
    auto [a, b, c] = x;
    REQUIRE(a == 1);
    REQUIRE(b == 2);
    REQUIRE(c == 3);
}

TEST_CASE("Vector3::generate_basis") {
    SECTION("simple case") {
        Vector3<double> x = {2, 0, 0};
        auto [v1, v2, v3] = x.generate_basis();
        
        REQUIRE(v1 == Vector3<double>{1, 0, 0});
        REQUIRE_THAT(v1.norm(), Catch::Matchers::WithinAbs(1.0, 1e-10));
        REQUIRE_THAT(v2.norm(), Catch::Matchers::WithinAbs(1.0, 1e-10));
        REQUIRE_THAT(v3.norm(), Catch::Matchers::WithinAbs(1.0, 1e-10));
        REQUIRE_THAT(v1.dot(v2), Catch::Matchers::WithinAbs(0.0, 1e-10));
        REQUIRE_THAT(v2.dot(v3), Catch::Matchers::WithinAbs(0.0, 1e-10));
        REQUIRE_THAT(v3.dot(v1), Catch::Matchers::WithinAbs(0.0, 1e-10));
    }
}
