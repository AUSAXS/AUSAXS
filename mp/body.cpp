// includes
#include <vector>
#include <string>
#include <fstream>

#include "Test.h"
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

// Test that the histograms are correct for proteins with only atoms (no waters)
void test_calc_distances_atoms() {
//*** TEST ATOMS ***//
    // the following just describes the eight corners of a cube centered at origo, with an additional atom at the very middle
    vector<Atom> a = {Atom(Vector3(-1, -1, -1), 1, "C", "C", 1), Atom(Vector3(-1, 1, -1), 1, "C", "C", 1),
                          Atom(Vector3(1, -1, -1), 1, "C", "C", 1), Atom(Vector3(1, 1, -1), 1, "C", "C", 1),
                          Atom(Vector3(-1, -1, 1), 1, "C", "C", 1), Atom(Vector3(-1, 1, 1), 1, "C", "C", 1),
                          Atom(Vector3(1, -1, 1), 1, "C", "C", 1), Atom(Vector3(1, 1, 1), 1, "C", "C", 1)};
    Body body(a, {});

    // set the weights to 1 so we can analytically determine the result
    for (auto& atom : body.protein_atoms) {
        atom.set_effective_charge(1);
    }
    body.updated_charge = true;

    // calculate the histogram
    shared_ptr<ScatteringHistogram> hist = body.get_histogram();
    const vector<double> d = hist->p_tot;

    // calculation: 8 identical points. 
    //      each point has:
    //          1 line  of length 0
    //          3 lines of length 2
    //          3 lines of length sqrt(2*2^2) = sqrt(8) = 2.82
    //          1 line  of length sqrt(3*2^2) = sqrt(12) = 3.46
    const vector<double> d_exp = {8, 0, 2*8*3, 8, 0, 0, 0, 0, 0, 0};
    for (size_t i = 0; i < d_exp.size(); i++) {
        if (d[i] != d_exp[i]) {
            IS_TRUE(false);
            cout << "Failed on index " << i << ". Values: " << d[i] << ", " << d_exp[i] << endl;
            // break;
        }
    }
}

// Test that the histograms are correct for proteins with only waters (no atoms)
void test_calc_distances_waters() {
    // the following just describes the eight corners of a cube centered at origo, with an additional atom at the very middle
    vector<Atom> a = {};
    vector<Hetatom> w = {Hetatom(Vector3(-1, -1, -1), 1, "C", "C", 1), Hetatom(Vector3(-1, 1, -1), 1, "C", "C", 1), 
                         Hetatom(Vector3(1, -1, -1), 1, "C", "C", 1),  Hetatom(Vector3(1, 1, -1), 1, "C", "C", 1), 
                         Hetatom(Vector3(-1, -1, 1), 1, "C", "C", 1),  Hetatom(Vector3(-1, 1, 1), 1, "C", "C", 1),
                         Hetatom(Vector3(1, -1, 1), 1, "C", "C", 1),   Hetatom(Vector3(1, 1, 1), 1, "C", "C", 1)};
    Body body(a, w);

    // set the weights to 1 so we can analytically determine the result
    for (auto& atom : body.hydration_atoms) {
        atom.set_effective_charge(1);
    }
    body.updated_charge = true;

    // calculate the histogram
    shared_ptr<ScatteringHistogram> hist = body.get_histogram();
    const vector<double> d = hist->p_tot;

    // calculation: 8 identical points. 
    //      each point has:
    //          1 line  of length 0
    //          3 lines of length 2
    //          3 lines of length sqrt(2*2^2) = sqrt(8) = 2.82
    //          1 line  of length sqrt(3*2^2) = sqrt(12) = 3.46
    const vector<double> d_exp = {8, 0, 2*8*3, 8, 0, 0, 0, 0, 0, 0};
    for (size_t i = 0; i < d_exp.size(); i++) {
        if (d[i] != d_exp[i]) {
            IS_TRUE(false);
            cout << "Failed on index " << i << ". Values: " << d[i] << ", " << d_exp[i] << endl;
            // break;
        }
    }
}

// Test that the histograms are correct for proteins with both atoms and waters
void test_calc_distances_both() {
    // the following just describes the eight corners of a cube centered at origo, with an additional atom at the very middle
    vector<Atom> a = {Atom(Vector3(-1, -1, -1), 1, "C", "C", 1), Atom(Vector3(-1, 1, -1), 1, "C", "C", 1),
                      Atom(Vector3(1, -1, -1), 1, "C", "C", 1), Atom(Vector3(1, 1, -1), 1, "C", "C", 1),
                      Atom(Vector3(-1, -1, 1), 1, "C", "C", 1), Atom(Vector3(-1, 1, 1), 1, "C", "C", 1)};
    vector<Hetatom> w = {Hetatom(Vector3(1, -1, 1), 1, "C", "C", 1),   Hetatom(Vector3(1, 1, 1), 1, "C", "C", 1)};
    Body body(a, w);

    // set the weights to 1 so we can analytically determine the result
    // waters
    for (auto& atom : body.hydration_atoms) {
        atom.set_effective_charge(1);
    }
    // atoms
    for (auto& atom : body.protein_atoms) {
        atom.set_effective_charge(1);
    }
    body.updated_charge = true;

    // calculate the histogram
    shared_ptr<ScatteringHistogram> hist = body.get_histogram();
    const vector<double> d = hist->p_tot;

    // calculation: 8 identical points. 
    //      each point has:
    //          1 line  of length 0
    //          3 lines of length 2
    //          3 lines of length sqrt(2*2^2) = sqrt(8) = 2.82
    //          1 line  of length sqrt(3*2^2) = sqrt(12) = 3.46
    const vector<double> d_exp = {8, 0, 2*8*3, 8, 0, 0, 0, 0, 0, 0};
    for (size_t i = 0; i < d_exp.size(); i++) {
        if (d[i] != d_exp[i]) {
            IS_TRUE(false);
            cout << "Failed on index " << i << ". Expected " << d_exp[i] << ", but received " << d[i] << endl;
            // break;
        }
    }
}

void test_translate() {
    Body body("temp.pdb");
    body.translate(Vector3({1, 1, 1}));
    IS_TRUE(body.protein_atoms[0].coords == Vector3({0, 0, 0}));
    IS_TRUE(body.protein_atoms[1].coords == Vector3({0, 2, 0}));
    IS_TRUE(body.protein_atoms[2].coords == Vector3({2, 0, 0}));
    IS_TRUE(body.protein_atoms[3].coords == Vector3({2, 2, 0}));
}

void test_rotate() {
    vector<Atom> a = {Atom(Vector3(1, 0, 0), 1, "C", "C", 1)};
    Body body(a, {});

    Vector3 axis = {0, 1, 0};
    body.rotate(axis, M_PI_2);
    IS_EQUAL(Vector3({0, 0, 1}), body.protein_atoms[0].coords); 
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

void test_update_charge() {
    Body body("temp.pdb");
    body.update_effective_charge();
}

int main(void) {
    cout << "Summary of Body testing:" << endl;
    create_test_file();
    test_get_cm();
    test_translate();
    test_rotate();
    test_calc_distances_atoms();
    test_calc_distances_waters();
    test_calc_distances_both();
    test_volume();
    test_get_mass();
    test_update_charge();
    remove("temp.pdb");

    if (passed_all) {
        cout << "\033[1;32m" << "All Body tests passed." << "\033[0m" << endl;
    } else {
        print_err("Some Body tests failed.");
    }
    return 0;
}