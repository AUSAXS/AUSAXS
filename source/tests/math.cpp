// includes
#include <vector>
#include <string>
#include <iostream>

#include "math/Matrix.h"
#include "math/Vector.h"
#include "math/Vector3.h"
#include "Tools.h"
#include "Test.h"

using std::cout, std::endl;

void test_vector3_basics() {
    Vector3 v = {1, 2, 3};
    Vector3 w = {4, 5, 6};

    IS_TRUE(v+w == Vector3({5, 7, 9}));
    IS_TRUE(v-w == Vector3({-3, -3, -3}));
    IS_TRUE(v.dot(w) == 4+10+18);
    IS_TRUE(v.norm() == 1+4+9);

    v += w; // v = (5, 7, 9)
    Vector3 a = w-v; // a = (-1, -2, -3)
    IS_TRUE(a == Vector3({-1, -2, -3}));
    v.x = 0;
    IS_TRUE(v == Vector3({0, 7, 9}));
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

int main(void) {
    cout << "Summary of math testing:" << endl;
    test_vector3_basics();
    test_matrix_basics();

    if (passed_all) {
        cout << "\033[1;32m" << "All math tests passed.           " << "\033[0m" << endl;
    } else {
        print_err("Some math tests failed.");
    }
    return 0;
}