#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>

#include <math/Matrix.h>
#include <math/Vector.h>
#include <math/Vector3.h>
#include <math/CubicSpline.h>
#include <math/LUPDecomposition.h>
#include <math/QRDecomposition.h>
#include <math/Statistics.h>
#include <dataset/SimpleDataset.h>
#include <plots/PlotDataset.h>

#include <vector>
#include <iostream>

using std::cout, std::endl;
using namespace ausaxs;

static double GenRandScalar() {
    return rand() % 100;
}

static Vector<double> GenRandVector(int m) {
    Vector<double> v(m);
    for (int i = 0; i < m; i++)
        v[i] = rand() % 100;
    return v;
}

static Matrix<double> GenRandMatrix(int n, int m) {
    Matrix<double> M(n, m);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            M[i][j] = rand() % 100;
    return M;
}

TEST_CASE("math: QRDecomposition") {
    Matrix<double> A = {{1, 2}, {3, 4}};
    QRDecomposition qr(A);
    auto I = matrix::identity(2);
    auto inv = qr.inverse();
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            REQUIRE_THAT((A*inv)[i][j], Catch::Matchers::WithinAbs(I[i][j], 1e-3));
        }
    }
    REQUIRE_THAT(qr.abs_determinant(), Catch::Matchers::WithinAbs(2, 1e-3));

    // randomized tests on 5x5 matrices
    srand(time(NULL)); // seed rng
    for (int i = 0; i < 10; i++) {
        A = GenRandMatrix(5, 5);
        Vector b = GenRandVector(5);
        QRDecomposition solver(A);
        Vector x = solver.solve(b);
        Vector Ax = A*x;
        auto id = matrix::identity(5);
        for (int j = 0; j < 5; j++) {
            REQUIRE_THAT(Ax[j], Catch::Matchers::WithinAbs(b[j], 1e-6));
            for (int k = 0; k < 5; k++) {
                REQUIRE_THAT((A*solver.inverse())[j][k], Catch::Matchers::WithinAbs(id[j][k], 1e-6));
            }
        }
    }
}

TEST_CASE("math: orthonormal_rotations") {
    for (int i = 0; i < 10; i++) {
        Vector3<double> angles = GenRandVector(3);
        Matrix R = matrix::rotation_matrix(angles.x(), angles.y(), angles.z());
        Matrix Ri = R.T();
        REQUIRE(R*Ri == matrix::identity(3));
    }

    for (int i = 0; i < 10; i++) {
        Vector3<double> axis = GenRandVector(3);
        double angle = GenRandScalar();
        Matrix R = matrix::rotation_matrix(axis, angle);
        Matrix Ri = R.T();
        REQUIRE(R*Ri == matrix::identity(3));
    }
}