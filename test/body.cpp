#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <vector>
#include <string>
#include <fstream>
#include <iostream>

#include <utility/Constants.h>
#include <hydrate/Grid.h>
#include <data/Protein.h>
#include <data/StateManager.h>
#include <data/BodySplitter.h>

using std::cout, std::endl, std::vector, std::shared_ptr;

TEST_CASE("translate", "[body]") {
    setting::protein::use_effective_charge = false;
    vector<Atom> a = {Atom(Vector3<double>(  -1, -1, -1), 1, "C", "C", 1), Atom(Vector3<double>(-1, 1, -1), 1, "C", "C", 1),
                        Atom(Vector3<double>( 1, -1, -1), 1, "C", "C", 1), Atom(Vector3<double>( 1, 1, -1), 1, "C", "C", 1),
                        Atom(Vector3<double>(-1, -1,  1), 1, "C", "C", 1), Atom(Vector3<double>(-1, 1,  1), 1, "C", "C", 1),
                        Atom(Vector3<double>( 1, -1,  1), 1, "C", "C", 1), Atom(Vector3<double>( 1, 1,  1), 1, "C", "C", 1)};
    Body body(a, {});

    SECTION("translate") {
        body.translate(Vector3<double>{1, 1, 1});
        CHECK(body.atoms(0).coords == Vector3<double>{0, 0, 0});
        CHECK(body.atoms(1).coords == Vector3<double>{0, 2, 0});
        CHECK(body.atoms(2).coords == Vector3<double>{2, 0, 0});
        CHECK(body.atoms(3).coords == Vector3<double>{2, 2, 0});
    }

    SECTION("trivial center") {
        body.translate(Vector3<double>{-1, -1, -1});
        body.center();
        CHECK(body.atoms(0).coords == Vector3<double>{-1, -1, -1});
        CHECK(body.atoms(1).coords == Vector3<double>{-1,  1, -1});
        CHECK(body.atoms(2).coords == Vector3<double>{ 1, -1, -1});
        CHECK(body.atoms(3).coords == Vector3<double>{ 1,  1, -1});
    }

    SECTION("non-trivial center") {
        a = {Atom(Vector3<double>(-1, -1, -1), 1, "C", "C", 1), Atom(Vector3<double>(-1, 1, -1), 1, "O", "O", 1),
             Atom(Vector3<double>( 1, -1, -1), 1, "C", "C", 1), Atom(Vector3<double>( 1, 1, -1), 1, "O", "O", 1)};
        body = Body(a, {});

        body.center();
        double shift = 0.142402;

        // make checks in a roundabout fashion for better failure messages
        auto res = Vector3<double>{-1, -1-shift, 0};
        if (body.atoms(0).coords.equals(res, 1e-3)) {CHECK(true);}
        else {CHECK(body.atoms(0).coords == res);}

        res = Vector3<double>{-1, 1-shift, 0};
        if (body.atoms(1).coords.equals(res, 1e-3)) {CHECK(true);}
        else {CHECK(body.atoms(1).coords == res);}

        res = Vector3<double>{1, -1-shift, 0};
        if (body.atoms(2).coords.equals(res, 1e-3)) {CHECK(true);}
        else {CHECK(body.atoms(2).coords == res);}

        res = Vector3<double>{1, 1-shift, 0};
        if (body.atoms(3).coords.equals(res, 1e-3)) {CHECK(true);}
        else {CHECK(body.atoms(3).coords == res);}
    }
}

TEST_CASE("body_rotate", "[body]") {
    SECTION("simple rotations") {
        vector<Atom> a = {Atom(Vector3<double>(1, 0, 0), 1, "C", "C", 1), 
                        Atom(Vector3<double>(0, 1, 0), 1, "C", "C", 1), 
                        Atom(Vector3<double>(0, 0, 1), 1, "C", "C", 1)};
        Body body(a, {});

        SECTION("axis") {
            Vector3<double> axis = {0, 1, 0};
            body.rotate(axis, M_PI_2);
            CHECK(Vector3<double>({0, 0, -1}) == body.atoms(0).coords); 
            CHECK(Vector3<double>({0, 1, 0}) == body.atoms(1).coords); 
            CHECK(Vector3<double>({1, 0, 0}) == body.atoms(2).coords); 

            axis = {1, 1, 1};
            body.rotate(axis, M_PI/4);
            CHECK(Vector3<double>({-0.5058793634, 0.3106172175, -0.8047378541}) == body.atoms(0).coords); 
            CHECK(Vector3<double>({-0.3106172175, 0.8047378541, 0.5058793634}) == body.atoms(1).coords); 
            CHECK(Vector3<double>({0.8047378541, 0.5058793634, -0.3106172175}) == body.atoms(2).coords); 
        }

        SECTION("euler") {
            body.rotate(0, M_PI_2, 0);
            CHECK(Vector3<double>({0, 0, -1}) == body.atoms(0).coords); 
            CHECK(Vector3<double>({0, 1, 0}) == body.atoms(1).coords); 
            CHECK(Vector3<double>({1, 0, 0}) == body.atoms(2).coords); 

            body.rotate(0.5612026, 0.3158423, 0.5612026);
            CHECK(Vector3<double>({-0.5058793634, 0.3106172175, -0.8047378541}).equals(body.atoms(0).coords, 1e-3));
            CHECK(Vector3<double>({-0.3106172175, 0.8047378541, 0.5058793634}).equals(body.atoms(1).coords, 1e-3));
            CHECK(Vector3<double>({0.8047378541, 0.5058793634, -0.3106172175}).equals(body.atoms(2).coords, 1e-3));
        }
    }

    SECTION("complex rotations") {
        vector<Atom> a = {Atom(Vector3<double>(0, 2, 1), 1, "C", "C", 1), 
                        Atom(Vector3<double>(5, 1, 3), 1, "C", "C", 1), 
                        Atom(Vector3<double>(6, 1, 4), 1, "C", "C", 1),
                        Atom(Vector3<double>(3, 7, 2), 1, "C", "C", 1)};
        Body body(a, {});

        Vector3<double> axis = {0.5, 2, 1};
        body.rotate(axis, 1.8);
        REQUIRE(Vector3<double>({0.5843819499, 1.6706126346, 1.3665837559}) == body.atoms(0).coords); 
        REQUIRE(Vector3<double>({1.8656722055, 4.7666664324, -2.9661689675}) == body.atoms(1).coords); 
        REQUIRE(Vector3<double>({2.6638285975, 5.6804357476, -3.692785794}) == body.atoms(2).coords); 
        REQUIRE(Vector3<double>({0.0886646879, 7.4409765368, 2.5737145825}) == body.atoms(3).coords); 
    }
}

TEST_CASE("body_get_cm", "[body]") {
    vector<Atom> a = {Atom(1, "C", "", "LYS", "", 1, "", Vector3<double>(-1, -1, -1), 1, 0, "C", "0"),  Atom(2, "C", "", "LYS", "", 1, "", Vector3<double>(-1, 1, -1), 1, 0, "C", "0"),
                        Atom(3, "C", "", "LYS", "", 1, "", Vector3<double>(1, -1, -1), 1, 0, "C", "0"), Atom(4, "C", "", "LYS", "", 1, "", Vector3<double>(1, 1, -1), 1, 0, "C", "0"),
                        Atom(5, "C", "", "LYS", "", 1, "", Vector3<double>(-1, -1, 1), 1, 0, "C", "0"), Atom(6, "C", "", "LYS", "", 1, "", Vector3<double>(-1, 1, 1), 1, 0, "C", "0"),
                        Atom(7, "C", "", "LYS", "", 1, "", Vector3<double>(1, -1, 1), 1, 0, "C", "0"),  Atom(8, "C", "", "LYS", "", 1, "", Vector3<double>(1, 1, 1), 1, 0, "C", "0")};
    Body body(a, {});
    Vector3<double> cm = body.get_cm();
    REQUIRE(cm == Vector3<double>({0, 0, 0}));
}

TEST_CASE("body_get_volume", "[body]") {
    vector<Atom> a = {Atom(1, "C", "", "LYS", "", 1, "", Vector3<double>(-1, -1, -1), 1, 0, "C", "0"),  Atom(2, "C", "", "LYS", "", 1, "", Vector3<double>(-1, 1, -1), 1, 0, "C", "0"),
                        Atom(3, "C", "", "LYS", "", 1, "", Vector3<double>(1, -1, -1), 1, 0, "C", "0"), Atom(4, "C", "", "LYS", "", 1, "", Vector3<double>(1, 1, -1), 1, 0, "C", "0"),
                        Atom(5, "C", "", "LYS", "", 1, "", Vector3<double>(-1, -1, 1), 1, 0, "C", "0"), Atom(6, "C", "", "LYS", "", 1, "", Vector3<double>(-1, 1, 1), 1, 0, "C", "0"),
                        Atom(7, "C", "", "LYS", "", 1, "", Vector3<double>(1, -1, 1), 1, 0, "C", "0"),  Atom(8, "C", "", "LYS", "", 1, "", Vector3<double>(1, 1, 1), 1, 0, "C", "0")};
    Body body(a, {});
    REQUIRE(body.get_volume_acids() == constants::volume::amino_acids.get("LYS"));
}

TEST_CASE("body_charge", "[body]") {
    vector<Atom> a = {Atom(Vector3<double>(  -1, -1, -1), 1, "C", "C", 1), Atom(Vector3<double>(-1, 1, -1), 1, "C", "C", 1),
                        Atom(Vector3<double>( 1, -1, -1), 1, "C", "C", 1), Atom(Vector3<double>( 1, 1, -1), 1, "C", "C", 1),
                        Atom(Vector3<double>(-1, -1,  1), 1, "C", "C", 1), Atom(Vector3<double>(-1, 1,  1), 1, "C", "C", 1),
                        Atom(Vector3<double>( 1, -1,  1), 1, "C", "C", 1), Atom(Vector3<double>( 1, 1,  1), 1, "C", "C", 1)};
    Body body(a, {});

    SECTION("effective charge") {
        double c0 = body.atoms(0).get_effective_charge();    
        double charge1 = body.total_effective_charge();
        body.update_effective_charge(1.5);
        double charge2 = body.total_effective_charge();
        CHECK(charge2 == charge1+12);
        CHECK(body.atoms(0).get_effective_charge() == c0+1.5);
    }

    SECTION("atomic charge") {
        CHECK(body.total_atomic_charge() == 8*6);
    }
}

TEST_CASE("body_mass", "[body]") {
    vector<Atom> a = {Atom(Vector3<double>(  -1, -1, -1), 1, "C", "C", 1), Atom(Vector3<double>(-1, 1, -1), 1, "C", "C", 1),
                        Atom(Vector3<double>( 1, -1, -1), 1, "C", "C", 1), Atom(Vector3<double>( 1, 1, -1), 1, "C", "C", 1),
                        Atom(Vector3<double>(-1, -1,  1), 1, "C", "C", 1), Atom(Vector3<double>(-1, 1,  1), 1, "C", "C", 1),
                        Atom(Vector3<double>( 1, -1,  1), 1, "C", "C", 1), Atom(Vector3<double>( 1, 1,  1), 1, "C", "C", 1)};
    Body body(a, {});

    CHECK_THAT(body.absolute_mass(), Catch::Matchers::WithinRel(8*constants::mass::atomic.get("C"), 1e-6));
    CHECK_THAT(body.molar_mass(), Catch::Matchers::WithinRel(8*constants::mass::atomic.get("C")*constants::Avogadro, 1e-6));
    utility::print_warning("Check definition of molar mass.");
}

TEST_CASE("body_equality", "[body]") {
    vector<Atom> a1 = {Atom(Vector3<double>(-1, -1, -1), 1, "C", "C", 1), Atom(Vector3<double>(-1, 1, -1), 1, "C", "C", 1)};
    vector<Atom> a2 = {Atom(Vector3<double>(-1, -1, -1), 1, "C", "C", 1), Atom(Vector3<double>(-1, 1, -1), 1, "C", "C", 1)};
    Body b1(a1);
    Body b2(a2);

    CHECK(!(b1 == b2)); // even though they have the same contents, body equality is defined exclusively by a uid
    Body b2c = b2;
    CHECK(b2 == b2c);
}

TEST_CASE("grid_add_remove_bodies", "[body]") {
    vector<Atom> a1 = {Atom(Vector3<double>(-1, -1, -1), 1, "C", "C", 1), Atom(Vector3<double>(-1, 1, -1), 1, "C", "C", 1)};
    vector<Atom> a2 = {Atom(Vector3<double>( 1, -1, -1), 1, "C", "C", 1), Atom(Vector3<double>( 1, 1, -1), 1, "C", "C", 1)};
    vector<Atom> a3 = {Atom(Vector3<double>(-1, -1,  1), 1, "C", "C", 1), Atom(Vector3<double>(-1, 1,  1), 1, "C", "C", 1)};
    vector<Atom> a4 = {Atom(Vector3<double>( 1, -1,  1), 1, "C", "C", 1), Atom(Vector3<double>( 1, 1,  1), 1, "C", "C", 1)};
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

TEST_CASE("split_body", "[body],[files]") {
    setting::general::verbose = false;
    vector<int> splits = {9, 99};
    Protein protein = BodySplitter::split("data/LAR1-2/LAR1-2.pdb", splits);

    // check sizes
    REQUIRE(protein.bodies.size() == 3);
    Body &b1 = protein.bodies[0], &b2 = protein.bodies[1], &b3 = protein.bodies[2];

    REQUIRE(b1.atoms().size() == 136);
    REQUIRE(b2.atoms().size() == 812-136);
    REQUIRE(b3.atoms().size() == 1606-812);

    // check start and end resseq
    CHECK(b1.atoms().back().resSeq == 8);
    CHECK(b2.atoms(0).resSeq == 9);
    CHECK(b2.atoms().back().resSeq == 98);
    CHECK(b3.atoms(0).resSeq == 99);
}

TEST_CASE("generate_sequential_constraints", "[body],[files]") {
    setting::general::verbose = false;
    vector<int> splits = {9, 99};
    Protein protein = BodySplitter::split("data/LAR1-2/LAR1-2.pdb", splits);
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

TEST_CASE("body_copy", "[body]") {
    vector<Atom> a1 = {Atom(Vector3<double>(-1, -1, -1), 1, "C", "C", 1), Atom(Vector3<double>(-1, 1, -1), 1, "C", "C", 1)};
    vector<Atom> a2 = {Atom(Vector3<double>( 1, -1, -1), 1, "C", "C", 1), Atom(Vector3<double>( 1, 1, -1), 1, "C", "C", 1)};
    vector<Atom> a3 = {Atom(Vector3<double>(-1, -1,  1), 1, "C", "C", 1), Atom(Vector3<double>(-1, 1,  1), 1, "C", "C", 1)};
    vector<Atom> a4 = {Atom(Vector3<double>( 1, -1,  1), 1, "C", "C", 1), Atom(Vector3<double>( 1, 1,  1), 1, "C", "C", 1)};
    Body b1(a1), b2(a2), b3(a3), b4(a4);

    SECTION("copy constructor") {
        Body b(b1);
        REQUIRE(b.atoms().size() == 2);
        REQUIRE(b.atoms(0) == a1[0]);
        REQUIRE(b.atoms(1) == a1[1]);

        // check that they are backed by separate files
        b.atoms(0) = a2[0];
        REQUIRE(b1.atoms(0) == a1[0]);
    }

    SECTION("move constructor") {
        Body b5 = std::move(b1);
        REQUIRE(b5.atoms().size() == 2);
        REQUIRE(b5.atoms(0) == a1[0]);
        REQUIRE(b5.atoms(1) == a1[1]);
    }

    SECTION("assignment") {
        b1 = b3;
        REQUIRE(b1.atoms().size() == 2);
        REQUIRE(b1.atoms(0) == a3[0]);
        REQUIRE(b1.atoms(1) == a3[1]);

        b1.atoms(0) = a1[0];
        REQUIRE(b3.atoms(0) == a3[0]);

        // assignment with temporary bodies
        b1 = Body();
        {
            Body b5(a1);
            b1 = b5;
        }
        REQUIRE(b1.atoms().size() == 2);
        REQUIRE(b1.atoms(0) == a1[0]);
        REQUIRE(b1.atoms(1) == a1[1]);
    }
}

TEST_CASE("body_io", "[body]") {
    vector<Atom> a = {Atom(1, "C"  , "", "LYS", "", 1, "", Vector3<double>(-1, -1, -1), 1, 0, "C", "0"),  Atom(2, "C", "", "LYS", "", 1, "", Vector3<double>(-1, 1, -1), 1, 0, "C", "0"),
                      Atom(3, "O"  , "", "LYS", "", 1, "", Vector3<double>(1, -1, -1), 1, 0, "O", "0"),   Atom(4, "C", "", "LYS", "", 1, "", Vector3<double>(1, 1, -1), 1, 0, "C", "0"),
                      Atom(5, "N"  , "", "LYS", "", 1, "", Vector3<double>(-1, -1, 1), 1, 0, "N", "0"),   Atom(6, "C", "", "LYS", "", 1, "", Vector3<double>(-1, 1, 1), 1, 0, "C", "0"),
                      Atom(7, "OXT", "", "LYS", "", 1, "", Vector3<double>(1, -1, 1), 1, 0, "O", "0"),    Atom(8, "C", "", "LYS", "", 1, "", Vector3<double>(1, 1, 1), 1, 0, "C", "0")};
    Body body(a, {});

    body.save("temp/body_io.pdb");
    Body body2("temp/body_io.pdb");

    CHECK(body.atoms().size() == body2.atoms().size());
    for (unsigned int i = 0; i < body.atoms().size(); i++) {
        CHECK(body.atoms(i).as_pdb() == body2.atoms(i).as_pdb());
    }

    remove("temp/body_io.pdb");
}