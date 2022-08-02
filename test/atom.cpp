#include <catch2/catch_test_macros.hpp>

#include "data/Atom.h"
#include "data/Hetatom.h"

TEST_CASE("Atom", "[atom]") {
    // "element", "resName", and "name" are used for some internal logic, and must have reasonable values. "" can also be used. 
//*** ATOMS ***//
    Atom a1 = Atom(15, "", "altLoc", "", "chainID", 3, "iCode", Vector3<double>({0, 1, 2}), 2.5, 3.5, "He", "2-");
    Atom a2 = Atom({3, 0, 5}, 2, "He", "", 3);

    CHECK(a1.serial == 15);
    CHECK(a1.name == "");
    CHECK(a1.altLoc == "altLoc");
    CHECK(a1.resName == "");
    CHECK(a1.chainID == "chainID");
    CHECK(a1.resSeq == 3);
    CHECK(a1.iCode == "iCode");
    CHECK(a1.coords == Vector3({0, 1, 2}));
    CHECK(a1.occupancy == 2.5);
    CHECK(a1.tempFactor == 3.5);
    CHECK(a1.element == "He");
    CHECK(a1.charge == "2-");

    CHECK(a2.serial == 3);
    CHECK(a2.name == "");
    CHECK(a2.altLoc == "");
    CHECK(a2.resName == "");
    CHECK(a2.chainID == "");
    CHECK(a2.resSeq == -1);
    CHECK(a2.iCode == "");
    CHECK(a2.coords == Vector3({3, 0, 5}));
    CHECK(a2.occupancy == 2);
    CHECK(a2.tempFactor == -1);
    CHECK(a2.element == "He");
    CHECK(a2.charge == "");

//*** HETATOMS ***//
    Hetatom w1 = Hetatom::create_new_water(Vector3<double>({1, 2, 3}));
    CHECK(w1.serial == -1);
    CHECK(w1.name == "O");
    CHECK(w1.altLoc == "");
    CHECK(w1.resName == "HOH");
    CHECK(w1.chainID == "");
    CHECK(w1.resSeq == -1);
    CHECK(w1.iCode == "");
    CHECK(w1.coords == Vector3({1, 2, 3}));
    CHECK(w1.occupancy == 1);
    CHECK(w1.tempFactor == 0);
    CHECK(w1.element == "O");
    CHECK(w1.charge == "");
}

TEST_CASE("operators", "[atom]") {
//*** ATOMS ***//
    Atom a1 = Atom({3, 0, 5}, 2, "He", "", 3);
    Atom a2 = a1;
    REQUIRE(a1 == a2);

    a2 = Atom({0, 4, 1}, 2, "He", "", 2);
    REQUIRE(a1 != a2);
    REQUIRE(a2 < a1);

//*** HETATOMS ***//
    Hetatom w1 = Hetatom({3, 0, 5}, 2, "He", "", 3);
    Hetatom w2 = w1;
    REQUIRE(w1 == w2);

    w2 = Atom({0, 4, 1}, 2, "He", "", 2);
    REQUIRE(w1 != w2);
    REQUIRE(w2 < w1);
}