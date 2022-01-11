// includes
#include <vector>
#include <string>
#include <fstream>

#include "data/Protein.h"
#include "hydrate/Grid.h"
#include "constants.h"
#include "data/StateManager.h"

#include "catch2/catch.hpp"

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
TEST_CASE("body_histogram", "[body]") {
    SECTION("atoms_only") {
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
                cout << "Failed on index " << i << ". Values: " << d[i] << ", " << d_exp[i] << endl;
                REQUIRE(false);
            }
        }
    }

    SECTION("waters_only") {
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
                cout << "Failed on index " << i << ". Values: " << d[i] << ", " << d_exp[i] << endl;
                REQUIRE(false);
            }
        }
    }

    SECTION("both_atoms_and_water") {
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
                cout << "Failed on index " << i << ". Expected " << d_exp[i] << ", but received " << d[i] << endl;
                REQUIRE(false);
            }
        }
    }
}

TEST_CASE("translate", "[body]") {
    vector<Atom> a = {Atom(Vector3(-1, -1, -1), 1, "C", "C", 1), Atom(Vector3(-1, 1, -1), 1, "C", "C", 1),
                        Atom(Vector3(1, -1, -1), 1, "C", "C", 1), Atom(Vector3(1, 1, -1), 1, "C", "C", 1),
                        Atom(Vector3(-1, -1, 1), 1, "C", "C", 1), Atom(Vector3(-1, 1, 1), 1, "C", "C", 1),
                        Atom(Vector3(1, -1, 1), 1, "C", "C", 1), Atom(Vector3(1, 1, 1), 1, "C", "C", 1)};
    Body body(a, {});

    body.translate(Vector3{1, 1, 1});
    REQUIRE(body.protein_atoms[0].coords == Vector3{0, 0, 0});
    REQUIRE(body.protein_atoms[1].coords == Vector3{0, 2, 0});
    REQUIRE(body.protein_atoms[2].coords == Vector3{2, 0, 0});
    REQUIRE(body.protein_atoms[3].coords == Vector3{2, 2, 0});
}

TEST_CASE("rotate", "[body]") {
    SECTION("simple rotations") {
        vector<Atom> a = {Atom(Vector3(1, 0, 0), 1, "C", "C", 1), 
                        Atom(Vector3(0, 1, 0), 1, "C", "C", 1), 
                        Atom(Vector3(0, 0, 1), 1, "C", "C", 1)};
        Body body(a, {});

        Vector3 axis = {0, 1, 0};
        body.rotate(axis, M_PI_2);
        REQUIRE(Vector3({0, 0, -1}) == body.protein_atoms[0].coords); 
        REQUIRE(Vector3({0, 1, 0}) == body.protein_atoms[1].coords); 
        REQUIRE(Vector3({1, 0, 0}) == body.protein_atoms[2].coords); 

        axis = {1, 1, 1};
        body.rotate(axis, M_PI/4);
        REQUIRE(Vector3({-0.5058793634, 0.3106172175, -0.8047378541}) == body.protein_atoms[0].coords); 
        REQUIRE(Vector3({-0.3106172175, 0.8047378541, 0.5058793634}) == body.protein_atoms[1].coords); 
        REQUIRE(Vector3({0.8047378541, 0.5058793634, -0.3106172175}) == body.protein_atoms[2].coords); 
    }

    SECTION("complex rotations") {
        vector<Atom> a = {Atom(Vector3(0, 2, 1), 1, "C", "C", 1), 
                        Atom(Vector3(5, 1, 3), 1, "C", "C", 1), 
                        Atom(Vector3(6, 1, 4), 1, "C", "C", 1),
                        Atom(Vector3(3, 7, 2), 1, "C", "C", 1)};
        Body body(a, {});

        Vector3 axis = {0.5, 2, 1};
        body.rotate(axis, 1.8);
        REQUIRE(Vector3({0.5843819499, 1.6706126346, 1.3665837559}) == body.protein_atoms[0].coords); 
        REQUIRE(Vector3({1.8656722055, 4.7666664324, -2.9661689675}) == body.protein_atoms[1].coords); 
        REQUIRE(Vector3({2.6638285975, 5.6804357476, -3.692785794}) == body.protein_atoms[2].coords); 
        REQUIRE(Vector3({0.0886646879, 7.4409765368, 2.5737145825}) == body.protein_atoms[3].coords); 
    }
}

TEST_CASE("body_get_mass", "[body]") {
    vector<Atom> a = {Atom(1, "C", "", "LYS", "", 1, "", Vector3(-1, -1, -1), 1, 0, "C", "0"),  Atom(2, "C", "", "LYS", "", 1, "", Vector3(-1, 1, -1), 1, 0, "C", "0"),
                        Atom(3, "C", "", "LYS", "", 1, "", Vector3(1, -1, -1), 1, 0, "C", "0"), Atom(4, "C", "", "LYS", "", 1, "", Vector3(1, 1, -1), 1, 0, "C", "0"),
                        Atom(5, "C", "", "LYS", "", 1, "", Vector3(-1, -1, 1), 1, 0, "C", "0"), Atom(6, "C", "", "LYS", "", 1, "", Vector3(-1, 1, 1), 1, 0, "C", "0"),
                        Atom(7, "C", "", "LYS", "", 1, "", Vector3(1, -1, 1), 1, 0, "C", "0"),  Atom(8, "C", "", "LYS", "", 1, "", Vector3(1, 1, 1), 1, 0, "C", "0")};
    Body body(a, {});

    REQUIRE(body.get_mass() == Approx(8*constants::mass::C));
}

TEST_CASE("body_get_cm", "[body]") {
    vector<Atom> a = {Atom(1, "C", "", "LYS", "", 1, "", Vector3(-1, -1, -1), 1, 0, "C", "0"),  Atom(2, "C", "", "LYS", "", 1, "", Vector3(-1, 1, -1), 1, 0, "C", "0"),
                        Atom(3, "C", "", "LYS", "", 1, "", Vector3(1, -1, -1), 1, 0, "C", "0"), Atom(4, "C", "", "LYS", "", 1, "", Vector3(1, 1, -1), 1, 0, "C", "0"),
                        Atom(5, "C", "", "LYS", "", 1, "", Vector3(-1, -1, 1), 1, 0, "C", "0"), Atom(6, "C", "", "LYS", "", 1, "", Vector3(-1, 1, 1), 1, 0, "C", "0"),
                        Atom(7, "C", "", "LYS", "", 1, "", Vector3(1, -1, 1), 1, 0, "C", "0"),  Atom(8, "C", "", "LYS", "", 1, "", Vector3(1, 1, 1), 1, 0, "C", "0")};
    Body body(a, {});
    Vector3 cm = body.get_cm();
    REQUIRE(cm == Vector3({0, 0, 0}));
}

TEST_CASE("body_get_volume", "[body]") {
    vector<Atom> a = {Atom(1, "C", "", "LYS", "", 1, "", Vector3(-1, -1, -1), 1, 0, "C", "0"),  Atom(2, "C", "", "LYS", "", 1, "", Vector3(-1, 1, -1), 1, 0, "C", "0"),
                        Atom(3, "C", "", "LYS", "", 1, "", Vector3(1, -1, -1), 1, 0, "C", "0"), Atom(4, "C", "", "LYS", "", 1, "", Vector3(1, 1, -1), 1, 0, "C", "0"),
                        Atom(5, "C", "", "LYS", "", 1, "", Vector3(-1, -1, 1), 1, 0, "C", "0"), Atom(6, "C", "", "LYS", "", 1, "", Vector3(-1, 1, 1), 1, 0, "C", "0"),
                        Atom(7, "C", "", "LYS", "", 1, "", Vector3(1, -1, 1), 1, 0, "C", "0"),  Atom(8, "C", "", "LYS", "", 1, "", Vector3(1, 1, 1), 1, 0, "C", "0")};
    Body body(a, {});
    REQUIRE(body.get_volume_acids() == constants::volume::lysine);
}

TEST_CASE("update_charge", "[broken],[body]") {
    vector<Atom> a = {Atom(Vector3(-1, -1, -1), 1, "C", "C", 1), Atom(Vector3(-1, 1, -1), 1, "C", "C", 1),
                        Atom(Vector3(1, -1, -1), 1, "C", "C", 1), Atom(Vector3(1, 1, -1), 1, "C", "C", 1),
                        Atom(Vector3(-1, -1, 1), 1, "C", "C", 1), Atom(Vector3(-1, 1, 1), 1, "C", "C", 1),
                        Atom(Vector3(1, -1, 1), 1, "C", "C", 1), Atom(Vector3(1, 1, 1), 1, "C", "C", 1)};
    Body body(a, {});
    body.update_effective_charge();
}