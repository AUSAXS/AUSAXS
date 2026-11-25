#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <math/Matrix.h>
#include <math/MatrixUtils.h>
#include <math/QRDecomposition.h>
#include <math/Vector.h>
#include <utility/Random.h>

#include <random>

using namespace ausaxs;

TEST_CASE("QRDecomposition::QRDecomposition") {
    SECTION("2x2 matrix") {
        Matrix<double> A = {{1, 2}, {3, 4}};
        QRDecomposition qr(A);
        
        CHECK(qr.Q.N == 2);
        CHECK(qr.Q.M == 2);
        CHECK(qr.R.N == 2);
        CHECK(qr.R.M == 2);
    }

    SECTION("3x3 matrix") {
        Matrix<double> A = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
        QRDecomposition qr(A);
        
        CHECK(qr.Q.N == 3);
        CHECK(qr.Q.M == 3);
        CHECK(qr.R.N == 3);
        CHECK(qr.R.M == 3);
    }
}

TEST_CASE("QRDecomposition::inverse") {
    SECTION("2x2 matrix") {
        Matrix<double> A = {{1, 2}, {3, 4}};
        QRDecomposition qr(A);
        auto I = matrix::identity(2);
        auto inv = qr.inverse();
        
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                REQUIRE_THAT((A*inv)[i][j], Catch::Matchers::WithinAbs(I[i][j], 1e-3));
            }
        }
    }

    SECTION("3x3 matrix") {
        Matrix<double> A = {{1, 0, 0}, {0, 2, 0}, {0, 0, 3}};
        QRDecomposition qr(A);
        auto I = matrix::identity(3);
        auto inv = qr.inverse();
        
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                REQUIRE_THAT((A*inv)[i][j], Catch::Matchers::WithinAbs(I[i][j], 1e-6));
            }
        }
    }
}

TEST_CASE("QRDecomposition::solve") {
    SECTION("simple 2x2 system") {
        Matrix<double> A = {{1, 2}, {3, 4}};
        Vector<double> b = {5, 11};
        
        QRDecomposition qr(A);
        Vector<double> x = qr.solve(b);
        Vector<double> Ax = A * x;
        
        for (int i = 0; i < 2; i++) {
            REQUIRE_THAT(Ax[i], Catch::Matchers::WithinAbs(b[i], 1e-6));
        }
    }

    SECTION("3x3 system") {
        Matrix<double> A = {{2, 0, 0}, {0, 3, 0}, {0, 0, 4}};
        Vector<double> b = {2, 6, 12};
        
        QRDecomposition qr(A);
        Vector<double> x = qr.solve(b);
        Vector<double> Ax = A * x;
        
        for (int i = 0; i < 3; i++) {
            REQUIRE_THAT(Ax[i], Catch::Matchers::WithinAbs(b[i], 1e-6));
        }
        
        REQUIRE_THAT(x[0], Catch::Matchers::WithinAbs(1, 1e-6));
        REQUIRE_THAT(x[1], Catch::Matchers::WithinAbs(2, 1e-6));
        REQUIRE_THAT(x[2], Catch::Matchers::WithinAbs(3, 1e-6));
    }
}

TEST_CASE("QRDecomposition::abs_determinant") {
    SECTION("2x2 matrix") {
        Matrix<double> A = {{1, 2}, {3, 4}};
        QRDecomposition qr(A);
        
        REQUIRE_THAT(qr.abs_determinant(), Catch::Matchers::WithinAbs(2, 1e-3));
    }

    SECTION("diagonal matrix") {
        Matrix<double> A = {{2, 0, 0}, {0, 3, 0}, {0, 0, 4}};
        QRDecomposition qr(A);
        
        REQUIRE_THAT(qr.abs_determinant(), Catch::Matchers::WithinAbs(24, 1e-6));
    }

    SECTION("identity matrix") {
        Matrix<double> A = matrix::identity(3);
        QRDecomposition qr(A);
        
        REQUIRE_THAT(qr.abs_determinant(), Catch::Matchers::WithinAbs(1, 1e-6));
    }
}

TEST_CASE("QRDecomposition: random 5x5 systems") {
    std::uniform_int_distribution<int> entry(0, 99);
    auto id = matrix::identity(5);
    for (int i = 0; i < 10; i++) {
        Matrix<double> A(5, 5);
        Vector<double> b(5);
        for (int j = 0; j < 5; j++) {
            b[j] = entry(random::generator());
            for (int k = 0; k < 5; k++) {
                A[j][k] = entry(random::generator());
            }
        }

        QRDecomposition solver(A);
        Vector<double> Ax = A*solver.solve(b);
        auto AAinv = A*solver.inverse();
        for (int j = 0; j < 5; j++) {
            REQUIRE_THAT(Ax[j], Catch::Matchers::WithinAbs(b[j], 1e-6));
            for (int k = 0; k < 5; k++) {
                REQUIRE_THAT(AAinv[j][k], Catch::Matchers::WithinAbs(id[j][k], 1e-6));
            }
        }
    }
}
