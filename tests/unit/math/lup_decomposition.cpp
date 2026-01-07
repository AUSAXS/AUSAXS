#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <math/LUPDecomposition.h>
#include <math/Matrix.h>
#include <math/MatrixUtils.h>

using namespace ausaxs;

TEST_CASE("LUPDecomposition::LUPDecomposition") {
    SECTION("2x2 matrix") {
        Matrix<double> A = {{1, 2}, {3, 4}};
        LUPDecomposition lup(A);
        lup.decompose();
        
        CHECK(lup.permutations >= 0);
    }

    SECTION("3x3 matrix") {
        Matrix<double> A = {{1, 2, 3}, {4, 5, 6}, {7, 8, 10}};
        LUPDecomposition lup(A);
        lup.decompose();
        
        CHECK(lup.permutations >= 0);
    }
}

TEST_CASE("LUPDecomposition::determinant") {
    SECTION("2x2 matrix") {
        Matrix<double> A = {{4, 1}, {2, 3}};
        LUPDecomposition lup(A);
        lup.decompose();
        
        double det = lup.determinant();
        CHECK_THAT(det, Catch::Matchers::WithinAbs(10, 0.6));
    }

    SECTION("identity matrix") {
        Matrix<double> A = matrix::identity(3);
        LUPDecomposition lup(A);
        lup.decompose();
        
        double det = lup.determinant();
        CHECK_THAT(det, Catch::Matchers::WithinAbs(1, 1e-6));
    }
}
