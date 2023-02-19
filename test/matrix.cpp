#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>

#include <vector>
#include <string>
#include <iostream>

#include <math/Matrix.h>
#include <math/Vector.h>
#include <math/Vector3.h>

TEST_CASE("basic operations") {
    Matrix<double> A({{1, 2}, {3, 4}});
    Matrix<double> B({{5, 6}, {7, 8}});
    Matrix<double> C({{2, 3}, {3, 4}});

    // addition
    REQUIRE(A+B == Matrix{{6, 8}, {10, 12}});
    REQUIRE(A+C == Matrix{{3, 5}, {6, 8}});
    REQUIRE(B+C == Matrix{{7, 9}, {10, 12}});

    // subtraction
    REQUIRE(A-B == Matrix{{-4, -4}, {-4, -4}});
    REQUIRE(A-C == Matrix{{-1, -1}, {0, 0}});
    REQUIRE(B-C == Matrix{{3, 3}, {4, 4}});

    REQUIRE(B-A == Matrix{{4, 4}, {4, 4}});
    REQUIRE(C-A == Matrix{{1, 1}, {0, 0}});
    REQUIRE(C-B == Matrix{{-3, -3}, {-4, -4}});

    // negation
    REQUIRE(-A == Matrix{{-1, -2}, {-3, -4}});
    REQUIRE(-B == Matrix{{-5, -6}, {-7, -8}});
    REQUIRE(-C == Matrix{{-2, -3}, {-3, -4}});

    // scalar multiplication
    REQUIRE(A*2 == Matrix{{2, 4}, {6, 8}});
    REQUIRE(B*3 == Matrix{{15, 18}, {21, 24}});
    REQUIRE(C*5 == Matrix{{10, 15}, {15, 20}});

    // scalar division
    REQUIRE(A/2 == Matrix{{1./2, 2./2}, {3./2, 4./2}});
    REQUIRE(B/3 == Matrix{{5./3, 6./3}, {7./3, 8./3}});
    REQUIRE(C/5 == Matrix{{2./5, 3./5}, {3./5, 4./5}});

    // matrix multiplication
    REQUIRE(A*B == Matrix{{19, 22}, {43, 50}});
    REQUIRE(A*C == Matrix{{8, 11}, {18, 25}});
    REQUIRE(B*C == Matrix{{28, 39}, {38, 53}});
}

TEST_CASE("assignment") {
    Matrix<double> A({{1, 2}, {3, 4}});
    Matrix<double> B({{5, 6}, {7, 8}});
    Matrix<double> C({{2, 3}, {3, 4}});

    A += C;
    B += C;
    REQUIRE(A == Matrix{{3, 5}, {6, 8}});
    REQUIRE(B == Matrix{{7, 9}, {10, 12}});

    A -= C;
    B -= C;
    REQUIRE(A == Matrix{{1, 2}, {3, 4}});
    REQUIRE(B == Matrix{{5, 6}, {7, 8}});

    A = C;
    B = C;
    C = Matrix<double>{{0, 0}, {0, 0}};
    REQUIRE(A == Matrix{{2, 3}, {3, 4}});
    REQUIRE(B == Matrix{{2, 3}, {3, 4}});

    A = B.copy();
    C = A.copy();
    B = matrix::identity(2);
    REQUIRE(A == Matrix{{2, 3}, {3, 4}});
    REQUIRE(C == A);
}

TEST_CASE("multiplication") {
    Matrix<double> A({{1, 2}, {3, 4}});
    Matrix<double> B({{5, 6}, {7, 8}});
    Matrix<double> C({{2, 3}, {3, 4}});

    // matrix multiplication
    B = Matrix<double>{{1, 2, 3}, {2, 3, 4}};
    REQUIRE(A*B == Matrix<double>{{5, 8, 11}, {11, 18, 25}});
    REQUIRE(C*B == Matrix<double>{{8, 13, 18}, {11, 18, 25}});

    // vector multiplication
    Vector<double> v = {1, 2};
    REQUIRE(A*v == Vector<double>{5, 11});
    REQUIRE(C*v == Vector<double>{8, 11});

    // transpose
    REQUIRE(A.T() == Matrix<double>{{1, 3}, {2, 4}});
    REQUIRE(B.T() == Matrix<double>{{1, 2}, {2, 3}, {3, 4}});
    REQUIRE(C.T() == Matrix<double>{{2, 3}, {3, 4}});
}

TEST_CASE("determinant") {
    Matrix<double> A({{1, 2}, {3, 4}});
    Matrix<double> B({{5, 6}, {7, 8}});
    Matrix<double> C({{2, 3}, {3, 4}});

    A = {{4, 1}, {2, 3}};
    B = {{-2, 3, -1}, {5, -1, 4}, {4, -8, 2}};
    C = {{5, -7, 2, 2}, {0, 3, 0, -4}, {-5, -8, 0, 3}, {0, 5, 0, -6}};
    REQUIRE_THAT(A.det(), Catch::Matchers::WithinAbs(10, 1e-3));
    REQUIRE_THAT(B.det(), Catch::Matchers::WithinAbs(-6, 1e-3));
    REQUIRE_THAT(C.det(), Catch::Matchers::WithinAbs(20, 1e-3));
}

TEST_CASE("rotations", "[math]") {
    Matrix<double> A({{1, 2}, {3, 4}});
    Matrix<double> B({{5, 6}, {7, 8}});
    Matrix<double> C({{2, 3}, {3, 4}});

    // check basic rotations
    Matrix R = matrix::rotation_matrix(M_PI/2, 0, 0);
    REQUIRE(R == Matrix{{1, 0, 0}, {0, 0, -1}, {0, 1, 0}});

    R = matrix::rotation_matrix(0, M_PI/2, 0);
    REQUIRE(R == Matrix{{0, 0, 1}, {0, 1, 0}, {-1, 0, 0}});

    R = matrix::rotation_matrix(0, 0, M_PI/2);
    REQUIRE(R == Matrix{{0, -1, 0}, {1, 0, 0}, {0, 0, 1}});
}