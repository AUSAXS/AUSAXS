// includes
#include <vector>
#include <string>
#include <fstream>

#include "tests/Test.h"
#include "data/Protein.h"
#include "hydrate/Grid.h"
#include "constants.h"
#include "data/StateManager.h"

void create_test_file() {
    std::ofstream file("temp.pdb");
    // the following just describes the eight corners of a cube centered at origo, with an additional atom at the very middle
    file << "ATOM      1  C   LYS A   1          -1      -1      -1  1.00 00.00           C \n"
         << "ATOM      2  C   LYS A   1          -1       1      -1  1.00 00.00           C \n"
         << "ATOM      3  C   LYS A   1           1      -1      -1  1.00 00.00           C \n"
         << "ATOM      4  C   LYS A   1           1       1      -1  1.00 00.00           C \n"

         << "ATOM      5  C   LYS A   1          -1      -1       1  1.00 00.00           C \n"
         << "ATOM      6  C   LYS A   1          -1       1       1  1.00 00.00           C \n"
         << "ATOM      7  C   LYS A   1           1      -1       1  1.00 00.00           C \n"
         << "ATOM      8  C   LYS A   1           1       1       1  1.00 00.00           C \n"

         << "ATOM      9  C   LYS A   1           0       0       0  1.00 00.00           C  ";
    file.close();
}

void test_translate() {
    Body body("temp.pdb");
    body.translate(Vector3({1, 1, 1}));
    IS_TRUE(body.protein_atoms[0].coords == Vector3({0, 0, 0}));
    IS_TRUE(body.protein_atoms[1].coords == Vector3({0, 2, 0}));
    IS_TRUE(body.protein_atoms[2].coords == Vector3({2, 0, 0}));
    IS_TRUE(body.protein_atoms[3].coords == Vector3({2, 2, 0}));
}

void test_get_mass() {
    Body body("temp.pdb");
    IS_TRUE(approx(body.get_mass(), 9*constants::mass::C));
}

void test_get_cm() {
    Body body("temp.pdb");
    Vector3 cm = body.get_cm();
    IS_TRUE(cm == Vector3({0, 0, 0}));
}

void test_volume() {
    Body body("temp.pdb");
    IS_TRUE(body.get_volume_acids() == constants::volume::lysine);
}

int main(void) {
    cout << "Summary of Body testing:" << endl;
    create_test_file();
    test_get_cm();
    test_translate();
    test_volume();
    test_get_mass();
    remove("temp.pdb");

    if (passed_all) {
        cout << "\033[1;32m" << "All Body tests passed." << "\033[0m" << endl;
    } else {
        print_err("Some Body tests failed.");
    }
    return 0;
}