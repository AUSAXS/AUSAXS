#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>

#include <math/Matrix.h>
#include <math/MatrixUtils.h>
#include <math/Vector.h>
#include <math/Vector3.h>

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

TEST_CASE("matrix::rotation_matrix is orthonormal") {
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
