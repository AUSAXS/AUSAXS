// includes
#include <vector>
#include <string>
#include <fstream>

#include "data/Protein.h"
#include "hydrate/Grid.h"
#include "constants.h"
#include "data/StateManager.h"
#include "data/BodySplitter.h"

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

TEST_CASE("grid_add_remove_bodies", "[body]") {
    vector<Atom> a1 = {Atom(Vector3(-1, -1, -1), 1, "C", "C", 1), Atom(Vector3(-1, 1, -1), 1, "C", "C", 1)};
    vector<Atom> a2 = {Atom(Vector3( 1, -1, -1), 1, "C", "C", 1), Atom(Vector3( 1, 1, -1), 1, "C", "C", 1)};
    vector<Atom> a3 = {Atom(Vector3(-1, -1,  1), 1, "C", "C", 1), Atom(Vector3(-1, 1,  1), 1, "C", "C", 1)};
    vector<Atom> a4 = {Atom(Vector3( 1, -1,  1), 1, "C", "C", 1), Atom(Vector3( 1, 1,  1), 1, "C", "C", 1)};
    Body b1(a1), b2(a2), b3(a3), b4(a4);
    vector<Body> bodies = {b1, b2, b3, b4};
    Protein protein(bodies);

    auto grid = protein.get_grid();
    REQUIRE(grid->a_members.size() == 8);
    grid->remove(&b2);
    REQUIRE(grid->a_members.size() == 6);
    grid->remove(&b1);
    REQUIRE(grid->a_members.size() == 4);
    grid->remove(&b3);
    REQUIRE(grid->a_members.size() == 2);

    auto remaining = grid->a_members;
    for (const auto& e : remaining) {
        REQUIRE((e == a4[0] || e == a4[1]));
    }

    // check volume
    REQUIRE(grid->volume != 0);
    grid->remove(&b4);
    REQUIRE(grid->volume == 0);
}

TEST_CASE("split_body", "[body]") {
    vector<int> splits = {9, 99};
    Protein protein = BodySplitter::split("data/LAR1-2.pdb", splits);

    // check sizes
    REQUIRE(protein.bodies.size() == 3);
    Body &b1 = protein.bodies[0], &b2 = protein.bodies[1], &b3 = protein.bodies[2];

    REQUIRE(b1.protein_atoms.size() == 136);
    REQUIRE(b2.protein_atoms.size() == 812-136);
    REQUIRE(b3.protein_atoms.size() == 1606-812);

    // check start and end resseq
    CHECK(b1.protein_atoms.back().resSeq == 8);
    CHECK(b2.protein_atoms[0].resSeq == 9);
    CHECK(b2.protein_atoms.back().resSeq == 98);
    CHECK(b3.protein_atoms[0].resSeq == 99);
}

TEST_CASE("generate_sequential_constraints", "[body]") {
    vector<int> splits = {9, 99};
    Protein protein = BodySplitter::split("data/LAR1-2.pdb", splits);
    vector<Constraint> constraints = BodySplitter::sequential_constraints(protein);

    REQUIRE(constraints.size() == 2);

    // check first constraint
    Constraint& c1 = constraints[0];
    REQUIRE(c1.atom1->name == "CA");
    REQUIRE(c1.atom1->serial == 131);
    REQUIRE(c1.atom2->name == "CA");
    REQUIRE(c1.atom2->serial == 138);

    // check second constraint
    Constraint& c2 = constraints[1];
    REQUIRE(c2.atom1->name == "CA");
    REQUIRE(c2.atom1->serial == 809);
    REQUIRE(c2.atom2->name == "CA");
    REQUIRE(c2.atom2->serial == 814);
}