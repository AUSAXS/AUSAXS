// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <math/Superposition.h>
#include <math/SymmetricEigen.h>
#include <math/MatrixUtils.h>
#include <math/Vector3.h>

#include <numbers>
#include <random>

using namespace ausaxs;
using namespace ausaxs::matrix;

namespace {
    std::vector<Vector3<double>> random_points(unsigned seed, int n = 15) {
        std::mt19937 gen(seed);
        std::uniform_real_distribution<double> dist(-8, 8);
        std::vector<Vector3<double>> v;
        v.reserve(n);
        for (int i = 0; i < n; ++i) {v.push_back({dist(gen), dist(gen), dist(gen)});}
        return v;
    }

    double matrix_diff(const Matrix<double>& A, const Matrix<double>& B) {
        double d = 0;
        for (unsigned int i = 0; i < 3; ++i) {
            for (unsigned int j = 0; j < 3; ++j) {d += std::abs(A(i, j) - B(i, j));}
        }
        return d;
    }
}

TEST_CASE("symmetric_eigen") {
    SECTION("diagonal matrix") {
        Matrix<double> A = {{3, 0, 0}, {0, 1, 0}, {0, 0, 2}};
        auto e = symmetric_eigen(A);
        CHECK_THAT(e.values[0], Catch::Matchers::WithinAbs(1, 1e-9));
        CHECK_THAT(e.values[1], Catch::Matchers::WithinAbs(2, 1e-9));
        CHECK_THAT(e.values[2], Catch::Matchers::WithinAbs(3, 1e-9));
    }

    SECTION("eigenpairs satisfy A v = lambda v") {
        Matrix<double> A = {{4, 1, -2}, {1, 2, 0}, {-2, 0, 3}};
        auto e = symmetric_eigen(A);
        for (unsigned int k = 0; k < 3; ++k) {
            Vector3<double> v{e.vectors[k][0], e.vectors[k][1], e.vectors[k][2]};
            Vector3<double> Av = A*v;
            Vector3<double> lv = e.values[k]*v;
            CHECK_THAT((Av - lv).magnitude(), Catch::Matchers::WithinAbs(0, 1e-8));
            CHECK_THAT(v.magnitude(), Catch::Matchers::WithinAbs(1, 1e-9)); // unit length
        }
    }

    SECTION("values are ascending") {
        Matrix<double> A = {{4, 1, -2}, {1, 2, 0}, {-2, 0, 3}};
        auto e = symmetric_eigen(A);
        CHECK(e.values[0] <= e.values[1]);
        CHECK(e.values[1] <= e.values[2]);
    }
}

TEST_CASE("optimal_rotation") {
    SECTION("recovers a known rotation") {
        auto R = matrix::rotation_matrix<double>(Vector3<double>{0.4, 1.0, -0.3}, 1.1);
        // optimal_rotation maximizes <R', R>, whose maximizer is R itself
        auto Ro = optimal_rotation(R);
        CHECK_THAT(matrix_diff(Ro, R), Catch::Matchers::WithinAbs(0, 1e-9));
    }

    SECTION("always returns a proper rotation") {
        // even for a matrix whose nearest orthogonal transform is a reflection, the result has det +1
        Matrix<double> M = {{-1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
        auto R = optimal_rotation(M);
        CHECK_THAT(R.det(), Catch::Matchers::WithinAbs(1, 1e-9));
    }
}

TEST_CASE("superpose") {
    SECTION("recovers an exact rigid transform") {
        auto axis = GENERATE(Vector3<double>{0,0,1}, Vector3<double>{1,0,0}, Vector3<double>{0.3,1,0.5});
        auto R = matrix::rotation_matrix<double>(axis, 0.9);
        Vector3<double> tr{2, -5, 3};

        auto from = random_points(1);
        std::vector<Vector3<double>> to;
        for (const auto& p : from) {to.push_back(R*p + tr);}

        auto res = superpose(from, to);
        CHECK_THAT(res.rmsd, Catch::Matchers::WithinAbs(0, 1e-9));
        CHECK_THAT(matrix_diff(res.rotation, R), Catch::Matchers::WithinAbs(0, 1e-8));
        CHECK_THAT((res.translation - tr).magnitude(), Catch::Matchers::WithinAbs(0, 1e-8));
    }

    SECTION("reports a positive residual for imperfect data") {
        auto from = random_points(2);
        auto to = random_points(3); // unrelated set
        auto res = superpose(from, to);
        CHECK(0 < res.rmsd);
    }

    SECTION("identity when the sets coincide") {
        auto from = random_points(4);
        auto res = superpose(from, from);
        CHECK_THAT(res.rmsd, Catch::Matchers::WithinAbs(0, 1e-9));
        CHECK_THAT(matrix_diff(res.rotation, Matrix<double>::identity(3)), Catch::Matchers::WithinAbs(0, 1e-8));
    }
}
