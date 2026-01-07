#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>

#include <math/Vector.h>

using namespace ausaxs;

TEST_CASE("Vector::Vector") {
    SECTION("empty") {
        Vector<double> x;
        REQUIRE(x.size() == 0);
    }

    SECTION("initializer_list") {
        Vector<double> x = {1, 2, 3, 4};
        REQUIRE(x.size() == 4);
        REQUIRE(x[0] == 1);
        REQUIRE(x[1] == 2);
        REQUIRE(x[2] == 3);
        REQUIRE(x[3] == 4);
    }

    SECTION("std::vector copy") {
        std::vector<double> v = {1, 2, 3};
        Vector<double> x(v);
        REQUIRE(x.size() == 3);
        REQUIRE(x[0] == 1);
        REQUIRE(x[1] == 2);
        REQUIRE(x[2] == 3);
    }

    SECTION("std::vector move") {
        std::vector<double> v = {1, 2, 3};
        Vector<double> x(std::move(v));
        REQUIRE(x.size() == 3);
        REQUIRE(x[0] == 1);
        REQUIRE(x[1] == 2);
        REQUIRE(x[2] == 3);
    }

    SECTION("size constructor") {
        Vector<double> x(5);
        REQUIRE(x.size() == 5);
        for (unsigned int i = 0; i < 5; ++i) {
            REQUIRE(x[i] == 0);
        }
    }

    SECTION("copy constructor") {
        Vector<double> x = {1, 2, 3};
        Vector<double> y(x);
        REQUIRE(y.size() == 3);
        REQUIRE(y == x);
    }

    SECTION("move constructor") {
        Vector<double> x = {1, 2, 3};
        Vector<double> y(std::move(x));
        REQUIRE(y.size() == 3);
        REQUIRE(y[0] == 1);
        REQUIRE(y[1] == 2);
        REQUIRE(y[2] == 3);
    }
}

TEST_CASE("Vector::operator=") {
    SECTION("copy assignment") {
        Vector<double> x = {1, 2, 3};
        Vector<double> y;
        y = x;
        REQUIRE(y == x);
    }

    SECTION("move assignment") {
        Vector<double> x = {1, 2, 3};
        Vector<double> y;
        y = std::move(x);
        REQUIRE(y.size() == 3);
        REQUIRE(y[0] == 1);
    }

    SECTION("initializer_list assignment") {
        Vector<double> x;
        x = {1, 2, 3};
        REQUIRE(x.size() == 3);
        REQUIRE(x[0] == 1);
    }
}

TEST_CASE("Vector::operator[]") {
    Vector<double> x = {1, 2, 3, 4};
    
    SECTION("const access") {
        const Vector<double>& cx = x;
        REQUIRE(cx[0] == 1);
        REQUIRE(cx[1] == 2);
        REQUIRE(cx[2] == 3);
        REQUIRE(cx[3] == 4);
    }

    SECTION("mutable access") {
        x[0] = 10;
        x[1] = 20;
        REQUIRE(x[0] == 10);
        REQUIRE(x[1] == 20);
    }
}

TEST_CASE("Vector::operator==") {
    Vector<double> x = {1, 2, 3, 4};
    Vector<double> y = {1, 2, 3, 4};
    Vector<double> z = {1, 2, 3, 5};

    REQUIRE(x == y);
    REQUIRE_FALSE(x == z);
}

TEST_CASE("Vector::operator!=") {
    Vector<double> x = {1, 2, 3, 4};
    Vector<double> y = {1, 2, 3, 4};
    Vector<double> z = {1, 2, 3, 5};

    REQUIRE_FALSE(x != y);
    REQUIRE(x != z);
}

TEST_CASE("Vector::operator+") {
    Vector<double> x = {1, 2, 3, 4};
    Vector<double> y = {2, 3, 4, 5};
    
    Vector<double> result = x + y;
    REQUIRE(result == Vector<double>{3, 5, 7, 9});
}

TEST_CASE("Vector::operator-") {
    Vector<double> x = {1, 2, 3, 4};
    Vector<double> y = {2, 3, 4, 5};
    
    Vector<double> result = x - y;
    REQUIRE(result == Vector<double>{-1, -1, -1, -1});
}

TEST_CASE("Vector::operator- (unary)") {
    Vector<double> x = {1, 2, 3, 4};
    Vector<double> result = -x;
    REQUIRE(result == Vector<double>{-1, -2, -3, -4});
}

TEST_CASE("Vector::operator* (vector)") {
    Vector<double> x = {1, 2, 3, 4};
    Vector<double> y = {2, 3, 4, 5};
    
    Vector<double> result = x * y;
    REQUIRE(result == Vector<double>{2, 6, 12, 20});
}

TEST_CASE("Vector::operator* (scalar)") {
    Vector<double> x = {1, 2, 3, 4};
    
    SECTION("right multiplication") {
        Vector<double> result = x * 2;
        REQUIRE(result == Vector<double>{2, 4, 6, 8});
    }

    SECTION("left multiplication") {
        Vector<double> result = 2 * x;
        REQUIRE(result == Vector<double>{2, 4, 6, 8});
    }
}

TEST_CASE("Vector::operator/") {
    Vector<double> x = {2, 4, 6, 8};
    Vector<double> result = x / 2;
    REQUIRE(result == Vector<double>{1, 2, 3, 4});
}

TEST_CASE("Vector::operator+=") {
    Vector<double> x = {1, 2, 3, 4};
    Vector<double> y = {2, 3, 4, 5};
    
    x += y;
    REQUIRE(x == Vector<double>{3, 5, 7, 9});
}

TEST_CASE("Vector::operator-=") {
    Vector<double> x = {1, 2, 3, 4};
    Vector<double> y = {2, 3, 4, 5};
    
    x -= y;
    REQUIRE(x == Vector<double>{-1, -1, -1, -1});
}

TEST_CASE("Vector::operator*=") {
    SECTION("scalar") {
        Vector<double> x = {1, 2, 3, 4};
        x *= 2;
        REQUIRE(x == Vector<double>{2, 4, 6, 8});
    }

    SECTION("vector") {
        Vector<double> x = {1, 2, 3, 4};
        Vector<double> y = {2, 3, 4, 5};
        x *= y;
        REQUIRE(x == Vector<double>{2, 6, 12, 20});
    }
}

TEST_CASE("Vector::operator/=") {
    Vector<double> x = {2, 4, 6, 8};
    x /= 2;
    REQUIRE(x == Vector<double>{1, 2, 3, 4});
}

TEST_CASE("Vector::dot") {
    Vector<double> x = {1, 2, 3, 4};
    Vector<double> y = {2, 3, 4, 5};
    
    double result = x.dot(y);
    REQUIRE(result == 2+6+12+20);
}

TEST_CASE("Vector::norm") {
    Vector<double> x = {1, 2, 3, 4};
    double result = x.norm();
    REQUIRE(result == sqrt(1+4+9+16));
}

TEST_CASE("Vector::norm2") {
    Vector<double> x = {1, 2, 3, 4};
    double result = x.norm2();
    REQUIRE(result == 1+4+9+16);
}

TEST_CASE("Vector::distance") {
    Vector<double> x = {1, 2, 3, 4};
    Vector<double> y = {2, 3, 4, 5};
    
    double result = x.distance(y);
    REQUIRE(result == sqrt(1+1+1+1));
}

TEST_CASE("Vector::distance2") {
    Vector<double> x = {1, 2, 3, 4};
    Vector<double> y = {2, 3, 4, 5};
    
    double result = x.distance2(y);
    REQUIRE(result == 1+1+1+1);
}

TEST_CASE("Vector::size") {
    Vector<double> x = {1, 2, 3, 4};
    REQUIRE(x.size() == 4);

    Vector<double> y;
    REQUIRE(y.size() == 0);
}

TEST_CASE("Vector::begin/end") {
    Vector<double> x = {1, 2, 3, 4};
    
    SECTION("const iterators") {
        const Vector<double>& cx = x;
        auto it = cx.begin();
        REQUIRE(*it == 1);
        ++it;
        REQUIRE(*it == 2);
        ++it;
        REQUIRE(*it == 3);
        ++it;
        REQUIRE(*it == 4);
        ++it;
        REQUIRE(it == cx.end());
    }

    SECTION("mutable iterators") {
        auto it = x.begin();
        *it = 10;
        REQUIRE(x[0] == 10);
    }

    SECTION("range-based for loop") {
        int count = 0;
        for (const auto& val : x) {
            count++;
            REQUIRE((val == 1 || val == 2 || val == 3 || val == 4));
        }
        REQUIRE(count == 4);
    }
}

TEST_CASE("Vector::conversion to std::vector") {
    Vector<double> x = {1, 2, 3, 4};
    std::vector<double> v = x;
    REQUIRE(v.size() == 4);
    REQUIRE(v[0] == 1);
    REQUIRE(v[1] == 2);
    REQUIRE(v[2] == 3);
    REQUIRE(v[3] == 4);
}
