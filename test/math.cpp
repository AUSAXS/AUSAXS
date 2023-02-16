#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>

#include <vector>
#include <string>
#include <iostream>
#include <random>

#include <math/Matrix.h>
#include <math/Vector.h>
#include <math/Vector3.h>
#include <math/Cramer2DSolver.h>
#include <math/GivensSolver.h>
#include <math/CubicSpline.h>
#include <math/LUPDecomposition.h>
#include <math/QRDecomposition.h>
#include <math/Statistics.h>
#include <dataset/SimpleDataset.h>
#include <plots/PlotDataset.h>

using std::cout, std::endl;

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

TEST_CASE("Cramer", "[math]") {
    Matrix<double> A = {{2, 3}, {3, -4}};
    Vector<double> b = {12, 1};
    Cramer2DSolver solver1(A);
    REQUIRE(solver1.solve(b) == Vector{3, 2});

    A = {{1, 2}, {4, 5}};
    b = {{3, 6}};
    Cramer2DSolver solver2(A);
    REQUIRE(solver2.solve(b) == Vector{-1, 2});

    A = {{2, -2}, {2, 2}};
    b = {{8, 2}};
    Cramer2DSolver solver3(A);
    REQUIRE(solver3.solve(b) == Vector{2.5, -1.5});

    // randomized tests on 2x2 matrices
    srand(time(NULL)); // seed rng
    for (int i = 0; i < 100; i++) {
        A = GenRandMatrix(2, 2);
        b = GenRandVector(2);
        Cramer2DSolver solver(A);
        Vector x = solver.solve(b);
        Vector Ax = A*x;
        REQUIRE(Ax == b);
    }
}

TEST_CASE("QRDecomposition", "[math],[broken]") {
    Matrix<double> A = {{1, 2}, {3, 4}};
    QRDecomposition qr(A);
    REQUIRE(A*qr.inverse() == matrix::identity(2));
    REQUIRE_THAT(qr.abs_determinant(), Catch::Matchers::WithinAbs(2, 1e-3));

    // randomized tests on 5x5 matrices
    srand(time(NULL)); // seed rng
    for (int i = 0; i < 10; i++) {
        A = GenRandMatrix(5, 5);
        Vector b = GenRandVector(5);
        QRDecomposition solver(A);
        Vector x = solver.solve(b);
        Vector Ax = A*x;
        REQUIRE(Ax == b);
        REQUIRE(solver.inverse()*A == matrix::identity(5));
    }
}

TEST_CASE("Givens", "[broken],[math]") {
    Matrix<double> A = {{2, 3}, {3, -4}};
    Vector<double> b = {12, 1};
    GivensSolver solver1(A);
    REQUIRE(solver1.solve(b) == Vector{3, 2});

    A = {{1, 2}, {4, 5}};
    b = {{3, 6}};
    GivensSolver solver2(A);
    REQUIRE(solver2.solve(b) == Vector{-1, 2});

    A = {{2, -2}, {2, 2}};
    b = {{8, 2}};
    GivensSolver solver3(A);
    REQUIRE(solver3.solve(b) == Vector{2.5, -1.5});

    A = {{2, 3, 4}, {5, -6, 7}, {8, 9, 10}};
    b = {{119, 80, 353}};
    GivensSolver solver4(A);
    REQUIRE(solver4.solve(b) == Vector{12, 13, 14});

    // randomized tests on 5x5 matrices
    // srand(time(NULL)); // seed rng
    // for (int i = 0; i < 10; i++) {
    //     A = GenRandMatrix(5, 5);
    //     b = GenRandVector(5);
    //     GivensSolver solver(A);
    //     Vector x = solver.solve(b);
    //     Vector Ax = A*x;
    //     Ax.print();
    //     b.print();
    //     std::cout << std::endl;
    //     IS_TRUE(Ax == b);
    // }
}

TEST_CASE("cubic_spline", "[manual],[math]") {
    double b = 2*M_PI; // the range is from 0 to this
    int len = 10;
    double step = b/len;
    Vector<double> x(len);
    Vector<double> y(len);
    for (int i = 0; i < len; i++) {
        x[i] = i*step;
        y[i] = sin(x[i]);
    }
    double steps = 4; // interpolates 4 steps between points

    CubicSpline csin(x, y);
    std::vector<double> newx, newy;

    for (size_t i = 0; i < x.size()-1; i++) {
        newx.push_back(x[i]);
        newy.push_back(y[i]);
        for (double z = x[i]; z < x[i+1]; z += (x[i+1] - x[i])/steps) {
            newx.push_back(z);
            newy.push_back(csin.spline(z));
        }
    }

    SimpleDataset original(x, y);
    SimpleDataset interpolated(newx, newy);
    original.add_plot_options(style::draw::points, {{"color", style::color::black}});
    interpolated.add_plot_options(style::draw::line, {{"color", style::color::red}});

    plots::PlotDataset plot(original);
    plot.plot(interpolated);
    plot.save("temp/cubicspline.png");
}

TEST_CASE("orthonormal_rotations", "[math]") {
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