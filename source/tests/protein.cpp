// includes
#include <vector>
#include <string>
#include <fstream>

#include "tests/Test.h"
#include "Protein.h"
#include "hydrate/Grid.h"
#include "data/properties.h"

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

void test_calc_distances() {
    std::ofstream file("temp_dist.pdb");
    // the following just describes the eight corners of a cube centered at origo, with an additional atom at the very middle
    file << "ATOM      1  C   LYS A   1          -1      -1      -1  1.00 00.00           C \n"
         << "ATOM      2  C   LYS A   1          -1       1      -1  1.00 00.00           C \n"
         << "ATOM      3  C   LYS A   1           1      -1      -1  1.00 00.00           C \n"
         << "ATOM      4  C   LYS A   1           1       1      -1  1.00 00.00           C \n"
         << "HETATM    5  HOH LYS A   1           1      -1      -1  1.00 00.00           O \n"
         << "HETATM    6  HOH LYS A   1           1       1      -1  1.00 00.00           O \n"
         << "HETATM    7  HOH LYS A   1           1       1      -1  1.00 00.00           O \n";
    file.close();

    Protein protein("temp_dist.pdb");
    remove("temp_dist.pdb");

    shared_ptr<Distances> d = protein.get_distances();
    IS_TRUE(d->d_pp.size() == 6);
    IS_TRUE(d->d_hh.size() == 3);
    IS_TRUE(d->d_hp.size() == 12);
}

void test_get_cm() {
    Protein protein("temp.pdb");
    Vector3 cm = protein.get_cm();
    IS_TRUE(abs(cm[0]) < 1e-9);
    IS_TRUE(abs(cm[1]) < 1e-9);
    IS_TRUE(abs(cm[2]) < 1e-9);
}

void test_volume() {
    Protein protein("temp.pdb");
    IS_TRUE(protein.get_volume_acids() == property::volume::lysine);
}

int main(void) {
    cout << "Summary of Protein testing:" << endl;
    create_test_file();
    test_get_cm();
    test_volume();
    remove("temp.pdb");

    if (passed_all) {
        cout << "\033[1;32m" << "All protein tests passed." << "\033[0m" << endl;
    } else {
        print_err("Some Protein tests failed.");
    }
    return 0;
}