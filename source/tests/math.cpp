// includes
#include <vector>
#include <string>
#include <iostream>

#include "math/VectorN.h"
#include "math/Vector3.h"
#include "Tools.h"
#include "Test.h"

using std::cout, std::endl;

void test_vector3() {
    Vector3 v = {1, 2, 3};
    Vector3 w = {4, 5, 6};

    IS_TRUE(v+w == Vector3({5, 7, 9}));
    IS_TRUE(v-w == Vector3({-3, -3, -3}));
    IS_TRUE(v.dot(w) == 4+10+18);
    IS_TRUE(v.norm() == 1+4+9);

    v += w; // v = (5, 7, 9)
    Vector3 a = w-v; // a = (-1, -2, -3)
    IS_TRUE(a == Vector3({-1, -2, -3}));
}

int main(void) {
    cout << "Summary of math testing:" << endl;
    test_vector3();

    if (passed_all) {
        cout << "\033[1;32m" << "All math tests passed.           " << "\033[0m" << endl;
    } else {
        print_err("Some math tests failed.");
    }
    return 0;
}