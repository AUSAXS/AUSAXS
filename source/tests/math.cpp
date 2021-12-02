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
#include "Tools.h"
#include "Test.h"

#include <TCanvas.h>
#include <TGraph.h>

using std::cout, std::endl;

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
    IS_TRUE(v.norm() == 1+4+9); // norm

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
    IS_TRUE(v*w == Vector({0, 4, 9, 8})); // vector multiplication
}

void test_matrix_basics() {
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
    C = A*B;
    IS_TRUE(C == Matrix({{5, 8, 11}, {11, 18, 25}}));

    // transpose
    IS_TRUE(C.T() == Matrix({{5, 11}, {8, 18}, {11, 25}}));
}

void test_matrix_advanced() {
    Matrix A = {{1, 2, 3}, {2, 3, 4}};
    Matrix B = A.copy();
    B[0][0] = 0;
    IS_TRUE(B != A); // copy is not by reference

    Vector v = {1, 2, 2};
    IS_TRUE(A*v == Vector({11, 16})); // vector multiplication
    IS_TRUE(B*v == Vector({10, 16})); // repeat

    {
        Matrix C = {{6, 5, 4}, {3, 2, 1}};
        B = C;
    }
    IS_TRUE(B == Matrix({{6, 5, 4}, {3, 2, 1}})); // assignment is not by reference
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
    Matrix x = solver4.solve(b);
    x.print();
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
    for (int i = 0; i < x.size()-1; i++) {
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
    test_vector3_basics();
    test_vector_advanced();
    test_matrix_basics();
    test_matrix_advanced();
    test_Cramer();
    test_Givens();
    test_cubic_spline();

    if (passed_all) {
        cout << "\033[1;32m" << "All math tests passed.           " << "\033[0m" << endl;
    } else {
        print_err("Some math tests failed.");
    }
    return 0;
}