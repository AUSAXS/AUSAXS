// includes
#include <vector>
#include <string>
#include <iostream>
#include <random>

#include "math/Matrix.h"
#include "math/Vector.h"
#include "math/Vector3.h"
#include "math/Cramer2DSolver.cpp"
#include "math/GivensSolver.cpp"
#include "math/CubicSpline.h"
#include "math/LUPDecomposition.h"
#include "math/QRDecomposition.h"
#include "Tools.h"
#include "Test.h"

#include <TCanvas.h>
#include <TGraph.h>

using std::cout, std::endl;

// test if two doubles are approximately equal
bool approx(double a, double b) {return a - 1e-9 < b && a + 1e-9 > b;}

class randomDouble {
public:
    randomDouble(double low, double high) : r(std::bind(std::uniform_real_distribution<>(low, high), std::default_random_engine())) {}
    double operator()(){return r(); }

private:
    std::function<double()> r;
};

Vector GenRandVector(int m) {
    randomDouble rand(0, 100);
    Vector v(m);
    for (int i = 0; i < m; i++)
        v[i] = rand();
    return v;
}

Matrix GenRandMatrix(int n, int m) {
    Matrix M(n, m);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            M[i][j] = rand();
    return M;
}

void test_vector3_basics() {
    Vector3 v = {1, 2, 3};
    Vector3 w = {4, 5, 6};

    IS_TRUE(v+w == Vector3({5, 7, 9})); // plus
    IS_TRUE(v-w == Vector3({-3, -3, -3})); // minus
    IS_TRUE(v.dot(w) == 4+10+18); // dot product
    IS_TRUE(v.norm() == sqrt(1+4+9)); // norm

    v += w; // v = (5, 7, 9)
    Vector3 a = w-v; // a = (-1, -2, -3)
    IS_TRUE(a == Vector3({-1, -2, -3}));
    v.x = 0;
    IS_TRUE(v == Vector3({0, 7, 9}));
}

void test_vector_advanced() {
    Vector v = {1, 2, 3, 4};
    Vector w = v.copy();
    v[0] = 0;
    IS_TRUE(w != v); // copy is not by reference

    {
        Vector x = {6, 5, 4, 3, 2, 1};
        w = x;
    }
    IS_TRUE(w == Vector({6, 5, 4, 3, 2, 1})); // assignment is not by reference

    w = {2, 2, 3, 3}; // list assignment
    IS_TRUE(v*w == Vector({0, 4, 9, 12})); // vector multiplication

    // iterator
    for (const auto& e : w) {
        IS_TRUE(e == 2 || e == 3);
    }
}

void test_matrix() {
    Matrix A({{1, 2}, {3, 4}});
    Matrix B({{5, 6}, {7, 8}});

    // addition/subtraction
    Matrix C = A + B;
    IS_TRUE(C == Matrix({{6, 8}, {10, 12}}));
    IS_TRUE(C - A == B);

    // multiplication
    C = A*B;
    IS_TRUE(C == Matrix({{19, 22}, {43, 50}}));

    // different shapes
    B = Matrix({{1, 2, 3}, {2, 3, 4}});
    C = A*B; // {{5, 8, 11}, {11, 18, 25}}
    IS_TRUE(C == Matrix({{5, 8, 11}, {11, 18, 25}}));

    // transpose
    IS_TRUE(C.T() == Matrix({{5, 11}, {8, 18}, {11, 25}}));

    // copy
    A = {{1, 2, 3}, {2, 3, 4}};
    B = A.copy();
    B[0][0] = 0;
    IS_TRUE(B != A); // copy is not by reference

    // vector multiplication
    Vector v = {1, 2, 2};
    IS_TRUE(A*v == Vector({11, 16}));
    IS_TRUE(B*v == Vector({10, 16}));

    {
        Matrix C = {{6, 5, 4}, {3, 2, 1}};
        B = C;
    }
    IS_TRUE(B == Matrix({{6, 5, 4}, {3, 2, 1}})); // assignment is not by reference

    // slice assignment
    B[0] = A[1];
    IS_TRUE(B == Matrix({{2, 3, 4}, {3, 2, 1}}));

    // determinants (based on LU decomposition)
    A = {{4, 1}, {2, 3}};
    B = {{-2, 3, -1}, {5, -1, 4}, {4, -8, 2}};
    Matrix C = {{5, -7, 2, 2}, {0, 3, 0, -4}, {-5, -8, 0, 3}, {0, 5, 0, -6}};
    IS_TRUE(approx(A.det(), 10));
    IS_TRUE(approx(B.det(), -6));
    IS_TRUE(approx(C.det(), 20));
}

void test_slices() {
    Matrix A = {{1, 1, 2, 2}, {3, 3, 2, 2}, {5, 5, 4, 4}};
    Matrix B;

    // row through operator[]
    IS_TRUE(A[0] == Vector({1, 1, 2, 2}));
    IS_TRUE(A[1] == Vector({3, 3, 2, 2}));
    IS_TRUE(A[2] == Vector({5, 5, 4, 4}));

    // explicit row
    IS_TRUE(A.row(0) == Vector({1, 1, 2, 2}));
    IS_TRUE(A.row(1) == Vector({3, 3, 2, 2}));
    IS_TRUE(A.row(2) == Vector({5, 5, 4, 4}));

    // col
    IS_TRUE(A.col(0) == Vector({1, 3, 5}));
    IS_TRUE(A.col(1) == Vector({1, 3, 5}));
    IS_TRUE(A.col(2) == Vector({2, 2, 4}));
    IS_TRUE(A.col(3) == Vector({2, 2, 4}));

    // assignment
    A.row(1) = {9, 1, 2, 3};
    A.row(2) = {6, 3, 1, 2};
    B = {{1, 1, 2, 2}, {9, 1, 2, 3}, {6, 3, 1, 2}};
    IS_TRUE(A == B);

    A.col(1) = {2, 5, 1};
    A.col(3) = {7, 1, 3};
    B = {{1, 2, 2, 7}, {9, 5, 2, 1}, {6, 1, 1, 3}};
    IS_TRUE(A == B);

    // minus-assignment
    A = B;
    A.row(0) -= A.row(1);
    A.row(1) -= A.row(2);
    IS_TRUE(A.row(0) == Vector({-8, -3, 0, 6}));
    IS_TRUE(A.row(1) == Vector({3, 4, 1, -2}));

    A = B;
    A.col(0) -= A.col(1);
    A.col(1) -= A.col(2);
    IS_TRUE(A.col(0) == Vector({-1, 4, 5}));
    IS_TRUE(A.col(1) == Vector({0, 3, 0}));

    // plus-assignment
    A = B;
    A.row(0) += A.row(1);
    A.row(1) += A.row(2);
    IS_TRUE(A.row(0) == Vector({10, 7, 4, 8}));
    IS_TRUE(A.row(1) == Vector({15, 6, 3, 4}));

    A = B;
    A.col(0) += A.col(1);
    A.col(1) += A.col(2);
    IS_TRUE(A.col(0) == Vector({3, 14, 7}));
    IS_TRUE(A.col(1) == Vector({4, 7, 2}));

    // vector cast
    A = {{1, 2, 2, 7}, {9, 5, 2, 1}, {6, 1, 1, 3}};
    Vector a = A.col(2);
    IS_TRUE(a.N == 3 && a.data.size() == 3);
    IS_TRUE(a == Vector({2, 2, 1}));
    IS_TRUE(A.col(2).operator Vector().N == 3 && A.col(2).operator Vector().data.size() == 3); // chain cast

    // dot with vector
    Vector b = {2, 3, 1, 5};
    IS_TRUE(A.row(0).dot(b) == (2+6+2+35));
    IS_TRUE(A.row(2).dot(b) == (12+3+1+15));

    b = {1, 4, 2};
    IS_TRUE(A.col(0).dot(b) == (1+36+12));
    IS_TRUE(A.col(2).dot(b) == (2+8+2));

    // dot with other slice
    IS_TRUE(A.col(0).dot(A.col(2)) == (2+18+6));
    IS_TRUE(A.row(0).dot(A.row(1)) == (9+10+4+7));

    // norm
    IS_TRUE(A.col(0).norm() == sqrt(1+81+36));
    IS_TRUE(A.row(0).norm() == sqrt(1+4+4+49));
}

void test_Cramer() {
    Matrix A = {{2, 3}, {3, -4}};
    Vector b = {12, 1};
    Cramer2DSolver solver1(A);
    IS_TRUE(solver1.solve(b) == Vector({3, 2}));

    A = {{1, 2}, {4, 5}};
    b = {{3, 6}};
    Cramer2DSolver solver2(A);
    IS_TRUE(solver2.solve(b) == Vector({-1, 2}));

    A = {{2, -2}, {2, 2}};
    b = {{8, 2}};
    Cramer2DSolver solver3(A);
    IS_TRUE(solver3.solve(b) == Vector({2.5, -1.5}));

    // randomized tests on 2x2 matrices
    for (int i = 0; i < 100; i++) {
        A = GenRandMatrix(2, 2);
        b = GenRandVector(2);
        Cramer2DSolver solver(A);
        Vector x = solver.solve(b);
        Vector Ax = A*x;
        IS_TRUE(Ax == b);
    }
}

void test_QRDecomposition() {
    Matrix A = {{1, 2}, {3, 4}};
    QRDecomposition qr(A);
    qr.Q.print();
    std::cout << "COMPARE: " << std::endl;
    std::cout << "	   0.316   0.949\n" << "          0.949  -0.316" << std::endl;
    IS_TRUE(A*qr.inverse() == Matrix({{1, 0}, {0, 1}}));
}

void test_Givens() {
    Matrix A = {{2, 3}, {3, -4}};
    Vector b = {12, 1};
    GivensSolver solver1(A);
    IS_TRUE(solver1.solve(b) == Vector({3, 2}));

    A = {{1, 2}, {4, 5}};
    b = {{3, 6}};
    GivensSolver solver2(A);
    IS_TRUE(solver2.solve(b) == Vector({-1, 2}));

    A = {{2, -2}, {2, 2}};
    b = {{8, 2}};
    GivensSolver solver3(A);
    IS_TRUE(solver3.solve(b) == Vector({2.5, -1.5}));

    A = {{2, 3, 4}, {5, -6, 7}, {8, 9, 10}};
    b = {{119, 80, 353}};
    GivensSolver solver4(A);
    // Matrix x = solver4.solve(b);
    // x.print();
    IS_TRUE(solver4.solve(b) == Vector({12, 13, 14}));

    // randomized tests on 5x5 matrices
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

void test_cubic_spline() {
    double b = 2*M_PI; // the range is from 0 to this
    int len = 10;
    double step = b/len;
    Vector x(len);
    Vector y(len);
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
    std::unique_ptr<TCanvas> c = std::make_unique<TCanvas>("c1", "canvas", 600, 600);
    std::unique_ptr<TGraph> g1 = std::make_unique<TGraph>(newx.size(), &newx[0], &newy[0]);
    std::unique_ptr<TGraph> g2 = std::make_unique<TGraph>(len, &x.data[0], &y.data[0]);
    g1->SetLineColor(kRed);
    g1->Draw("AC");
    g2->Draw("SAME *");

    // setup the canvas and save the plot
    string path = "temp/cubicspline.pdf";
    c->SetRightMargin(0.15);
    c->SetLeftMargin(0.15);
    c->SaveAs(path.c_str());
}

int main(void) {
    cout << "Summary of math testing:" << endl;
    test_vector3();
    test_vector();
    test_matrix();
    test_slices();
    test_QRDecomposition();
    // test_Cramer();
    // test_Givens();
    // test_cubic_spline();

    if (passed_all) {
        cout << "\033[1;32m" << "All math tests passed.           " << "\033[0m" << endl;
    } else {
        print_err("Some math tests failed.");
    }
    return 0;
}