#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>

#include <math/Matrix.h>
#include <math/Vector.h>
#include <math/MatrixUtils.h>

using namespace ausaxs;

TEST_CASE("Matrix::Matrix") {
    SECTION("default") {
        Matrix<double> A;
        CHECK(A.N == 0);
        CHECK(A.M == 0);
        CHECK(A.data.empty());
    }

    SECTION("copy constructor") {
        Matrix<double> A({{1, 2}, {3, 4}});
        Matrix<double> B(A);
        CHECK(B == A);
    }

    SECTION("move constructor") {
        Matrix<double> A({{1, 2}, {3, 4}});
        Matrix<double> B(std::move(A));
        CHECK(B.N == 2);
        CHECK(B.M == 2);
        CHECK(B.row(0) == std::vector{1, 2});
        CHECK(B.row(1) == std::vector{3, 4});
    }

    SECTION("initializer_list") {
        Matrix<double> A({{1, 2}, {3, 4}, {5, 6}});
        CHECK(A.N == 3);
        CHECK(A.M == 2);
        CHECK(A.row(0) == std::vector{1, 2});
        CHECK(A.row(1) == std::vector{3, 4});
        CHECK(A.row(2) == std::vector{5, 6});
    }

    SECTION("vector<vector>") {
        std::vector<std::vector<double>> data = {{1, 2}, {3, 4}, {5, 6}};
        Matrix<double> A(data);
        CHECK(A.N == 3);
        CHECK(A.M == 2);
        CHECK(A.row(0) == std::vector{1, 2});
        CHECK(A.row(1) == std::vector{3, 4});
        CHECK(A.row(2) == std::vector{5, 6});
    }

    SECTION("from Vector") {
        Vector<double> v = {1, 2, 3};
        Matrix<double> A(v);
        CHECK(A.N == 3);
        CHECK(A.M == 1);
        CHECK(A.row(0) == std::vector{1});
        CHECK(A.row(1) == std::vector{2});
        CHECK(A.row(2) == std::vector{3});
    }

    SECTION("size constructor") {
        Matrix<double> A(3, 2);
        CHECK(A.N == 3);
        CHECK(A.M == 2);
        CHECK(A.row(0) == std::vector{0, 0});
        CHECK(A.row(1) == std::vector{0, 0});
        CHECK(A.row(2) == std::vector{0, 0});
    }
}

TEST_CASE("Matrix::operator=") {
    SECTION("copy assignment") {
        Matrix<double> A({{1, 2}, {3, 4}});
        Matrix<double> B;
        B = A;
        REQUIRE(B == A);
    }

    SECTION("move assignment") {
        Matrix<double> A({{1, 2}, {3, 4}});
        Matrix<double> B;
        B = std::move(A);
        REQUIRE(B.N == 2);
        REQUIRE(B.M == 2);
    }
}

TEST_CASE("Matrix::operator[]") {
    Matrix<double> A({{1, 2}, {3, 4}});
    
    SECTION("const access") {
        const Matrix<double>& cA = A;
        CHECK(cA[0] == std::vector{1, 2});
        CHECK(cA[1] == std::vector{3, 4});
    }

    SECTION("mutable access") {
        A[0][0] = 10;
        CHECK(A[0][0] == 10);
    }
}

TEST_CASE("Matrix::operator==") {
    Matrix<double> A({{1, 2}, {3, 4}});
    Matrix<double> B({{1, 2}, {3, 4}});
    Matrix<double> C({{1, 2}, {3, 5}});

    REQUIRE(A == B);
    REQUIRE_FALSE(A == C);
}

TEST_CASE("Matrix::operator!=") {
    Matrix<double> A({{1, 2}, {3, 4}});
    Matrix<double> B({{1, 2}, {3, 4}});
    Matrix<double> C({{1, 2}, {3, 5}});

    REQUIRE_FALSE(A != B);
    REQUIRE(A != C);
}

TEST_CASE("Matrix::operator+") {
    Matrix<double> A({{1, 2}, {3, 4}});
    Matrix<double> B({{5, 6}, {7, 8}});
    
    Matrix<double> result = A + B;
    REQUIRE(result == Matrix<double>{{6, 8}, {10, 12}});
}

TEST_CASE("Matrix::operator-") {
    SECTION("binary") {
        Matrix<double> A({{1, 2}, {3, 4}});
        Matrix<double> B({{5, 6}, {7, 8}});
        
        Matrix<double> result = A - B;
        REQUIRE(result == Matrix<double>{{-4, -4}, {-4, -4}});
    }

    SECTION("unary") {
        Matrix<double> A({{1, 2}, {3, 4}});
        Matrix<double> result = -A;
        REQUIRE(result == Matrix<double>{{-1, -2}, {-3, -4}});
    }
}

TEST_CASE("Matrix::operator* (scalar)") {
    Matrix<double> A({{1, 2}, {3, 4}});
    
    SECTION("right multiplication") {
        Matrix<double> result = A * 2;
        REQUIRE(result == Matrix<double>{{2, 4}, {6, 8}});
    }

    SECTION("left multiplication") {
        Matrix<double> result = 2 * A;
        REQUIRE(result == Matrix<double>{{2, 4}, {6, 8}});
    }
}

TEST_CASE("Matrix::operator/ (scalar)") {
    Matrix<double> A({{2, 4}, {6, 8}});
    Matrix<double> result = A / 2;
    REQUIRE(result == Matrix<double>{{1, 2}, {3, 4}});
}

TEST_CASE("Matrix::operator* (matrix)") {
    Matrix<double> A({{1, 2}, {3, 4}});
    Matrix<double> B({{5, 6}, {7, 8}});
    
    Matrix<double> result = A * B;
    REQUIRE(result == Matrix<double>{{19, 22}, {43, 50}});
}

TEST_CASE("Matrix::operator* (vector)") {
    Matrix<double> A({{1, 2}, {3, 4}});
    Vector<double> v = {1, 2};
    
    Vector<double> result = A * v;
    REQUIRE(result == Vector<double>{5, 11});
}

TEST_CASE("Matrix::operator+=") {
    Matrix<double> A({{1, 2}, {3, 4}});
    Matrix<double> B({{5, 6}, {7, 8}});
    
    A += B;
    REQUIRE(A == Matrix<double>{{6, 8}, {10, 12}});
}

TEST_CASE("Matrix::operator-=") {
    Matrix<double> A({{1, 2}, {3, 4}});
    Matrix<double> B({{5, 6}, {7, 8}});
    
    A -= B;
    REQUIRE(A == Matrix<double>{{-4, -4}, {-4, -4}});
}

TEST_CASE("Matrix::operator*=") {
    Matrix<double> A({{1, 2}, {3, 4}});
    A *= 2;
    REQUIRE(A == Matrix<double>{{2, 4}, {6, 8}});
}

TEST_CASE("Matrix::operator/=") {
    Matrix<double> A({{2, 4}, {6, 8}});
    A /= 2;
    REQUIRE(A == Matrix<double>{{1, 2}, {3, 4}});
}

TEST_CASE("Matrix::row") {
    Matrix<double> A({{1, 2, 3}, {4, 5, 6}, {7, 8, 9}});
    
    SECTION("const access") {
        const Matrix<double>& cA = A;
        CHECK(cA.row(0) == std::vector{1, 2, 3});
        CHECK(cA.row(1) == std::vector{4, 5, 6});
        CHECK(cA.row(2) == std::vector{7, 8, 9});
    }

    SECTION("mutable access") {
        auto r = A.row(0);
        r[0] = 10;
        CHECK(A[0][0] == 10);
    }
}

TEST_CASE("Matrix::col") {
    Matrix<double> A({{1, 2, 3}, {4, 5, 6}, {7, 8, 9}});
    
    CHECK(A.col(0) == std::vector{1, 4, 7});
    CHECK(A.col(1) == std::vector{2, 5, 8});
    CHECK(A.col(2) == std::vector{3, 6, 9});
}

TEST_CASE("Matrix::T") {
    Matrix<double> A({{1, 2}, {3, 4}});
    Matrix<double> result = A.T();
    REQUIRE(result == Matrix<double>{{1, 3}, {2, 4}});
}

TEST_CASE("Matrix::det") {
    SECTION("2x2") {
        Matrix<double> A({{4, 1}, {2, 3}});
        REQUIRE_THAT(A.det(), Catch::Matchers::WithinAbs(10, 1e-3));
    }

    SECTION("3x3") {
        Matrix<double> B({{-2, 3, -1}, {5, -1, 4}, {4, -8, 2}});
        REQUIRE_THAT(B.det(), Catch::Matchers::WithinAbs(-6, 1e-3));
    }

    SECTION("4x4") {
        Matrix<double> C({{5, -7, 2, 2}, {0, 3, 0, -4}, {-5, -8, 0, 3}, {0, 5, 0, -6}});
        REQUIRE_THAT(C.det(), Catch::Matchers::WithinAbs(20, 1e-3));
    }
}

TEST_CASE("Matrix::identity") {
    SECTION("2x2") {
        Matrix<double> I = matrix::identity(2);
        REQUIRE(I == Matrix<double>{{1, 0}, {0, 1}});
    }

    SECTION("3x3") {
        Matrix<double> I = matrix::identity(3);
        REQUIRE(I == Matrix<double>{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}});
    }
}

TEST_CASE("Matrix::copy") {
    Matrix<double> A({{1, 2}, {3, 4}});
    Matrix<double> B = A.copy();
    
    REQUIRE(B == A);
    B[0][0] = 10;
    REQUIRE(A[0][0] == 1);
}

TEST_CASE("Matrix::push_back") {
    Matrix<double> A({{1, 2}, {3, 4}});
    A.push_back({5, 6});
    
    REQUIRE(A.N == 3);
    REQUIRE(A.M == 2);
    REQUIRE(A.row(2) == std::vector{5, 6});
}

TEST_CASE("Matrix::extend") {
    Matrix<double> A({{1, 2}, {3, 4}});
    A.extend(2);
    
    REQUIRE(A.N == 4);
    REQUIRE(A.M == 2);
    REQUIRE(A.row(0) == std::vector{1, 2});
    REQUIRE(A.row(1) == std::vector{3, 4});
    REQUIRE(A.row(2) == std::vector{0, 0});
    REQUIRE(A.row(3) == std::vector{0, 0});
}

TEST_CASE("Matrix::resize") {
    Matrix<double> A({{1, 2}, {3, 4}});
    A.resize(3, 3);
    
    REQUIRE(A.N == 3);
    REQUIRE(A.M == 3);
}

TEST_CASE("Matrix::begin/end") {
    Matrix<double> A({{1, 2}, {3, 4}});
    
    SECTION("const iterators") {
        const Matrix<double>& cA = A;
        auto it = cA.begin();
        REQUIRE(*it == 1);
        ++it;
        REQUIRE(*it == 2);
        ++it;
        REQUIRE(*it == 3);
        ++it;
        REQUIRE(*it == 4);
        ++it;
        REQUIRE(it == cA.end());
    }

    SECTION("mutable iterators") {
        auto it = A.begin();
        *it = 10;
        REQUIRE(A[0][0] == 10);
    }
}
