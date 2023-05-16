#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>

#include <data/Atom.h>
#include <data/Water.h>
#include <settings/All.h>
#include <utility/Constants.h>

TEST_CASE("constructors") {
    settings::protein::use_effective_charge = false;

    SECTION("atom") {
        // "element", "resName", and "name" are used for some internal logic, and must have reasonable values. "" can also be used. 
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
        CHECK(a1.get_type() == RecordType::ATOM);
        CHECK(a1.is_water() == false);

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
        CHECK(a2.get_type() == RecordType::ATOM);
        CHECK(a2.is_water() == false);
    }

    SECTION("water") {
        Water w1 = Water::create_new_water(Vector3<double>({1, 2, 3}));
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
        CHECK(w1.get_type() == RecordType::WATER);
        CHECK(w1.is_water() == true);
    }
}

TEST_CASE("setters_getters") {
    Atom a1;
    a1.set_x(2);
    a1.set_y(3);
    a1.set_z(4);
    a1.set_occupancy(2.5);
    a1.set_tempFactor(3.5);
    a1.set_altLoc("altLoc");
    a1.set_serial(15);
    a1.set_resSeq(3);
    a1.set_iCode("iCode");
    a1.set_chainID("chainID");
    a1.set_element("He");
    a1.set_charge("2-");
    a1.set_resName("resName");
    a1.set_name("name");

    CHECK(a1.get_serial() == 15);
    CHECK(a1.get_name() == "name");
    CHECK(a1.get_altLoc() == "altLoc");
    CHECK(a1.get_resName() == "resName");
    CHECK(a1.get_chainID() == "chainID");
    CHECK(a1.get_resSeq() == 3);
    CHECK(a1.get_iCode() == "iCode");
    CHECK(a1.get_coordinates() == Vector3({2, 3, 4}));
    CHECK(a1.get_occupancy() == 2.5);
    CHECK(a1.get_tempFactor() == 3.5);
    CHECK(a1.get_element() == "He");
    CHECK(a1.get_charge() == "2-");

    CHECK(a1.Z() == 2);
    CHECK(a1.get_mass() == constants::mass::atomic.get("He"));

    const auto a2 = a1;
    CHECK(a2.get_coordinates() == Vector3({2, 3, 4}));
}

TEST_CASE("use_effective_charge") {
    settings::protein::use_effective_charge = true;
    Atom a(15, "O", "altLoc", "LYS", "chainID", 3, "iCode", Vector3<double>({0, 1, 2}), 2.5, 3.5, "O", "0+");

    auto res = constants::charge::atomic.get("O") + constants::hydrogen_atoms::residues.get("LYS").get("O", "O");
    CHECK(a.get_effective_charge() == res);
    a.add_effective_charge(1.5);
    CHECK(a.get_effective_charge() == res+1.5);

    CHECK(a.get_mass() == constants::mass::atomic.get("O") + constants::hydrogen_atoms::residues.get("LYS").get("O", "O"));
    settings::protein::use_effective_charge = false;
    CHECK(a.get_mass() == constants::mass::atomic.get("O"));
}

TEST_CASE("pdb") {
    std::string s = "ATOM      1  CB  ARG A 129         2.1     3.2     4.3  0.50 42.04           C ";
    Atom a; 
    a.parse_pdb(s);
    CHECK(a.serial == 1);
    CHECK(a.name == "CB");
    CHECK(a.altLoc == " "); // spaces are only removed from number strings
    CHECK(a.resName == "ARG");
    CHECK(a.chainID == "A");
    CHECK(a.resSeq == 129);
    CHECK(a.iCode == " "); // same
    CHECK(a.coords == Vector3({2.1, 3.2, 4.3}));
    CHECK(a.occupancy == 0.5);
    CHECK(a.tempFactor == 42.04);
    CHECK(a.element == "C");
    CHECK(a.charge == "  "); // same
    CHECK(a.get_recName() == "ATOM  ");

    Atom b;
    b.parse_pdb(a.as_pdb());
    CHECK(a.equals_content(b));
}

TEST_CASE("coords") {
    Atom a1({0, 0, 0}, 1, "O", "", 1);
    Atom a2({1, 0, 0}, 1, "O", "", 1);
    Atom a3({0, 1, 0}, 1, "O", "", 1);
    Atom a4({0, 0, 1}, 1, "O", "", 1);
    Atom a5({1, 1, 1}, 1, "O", "", 1);

    CHECK(a1.distance(a2) == 1);
    CHECK(a1.distance(a3) == 1);
    CHECK(a1.distance(a4) == 1);
    CHECK_THAT(a1.distance(a5), Catch::Matchers::WithinAbs(std::sqrt(3), 1e-6));

    a1.translate({1, 1, 1});
    CHECK_THAT(a1.distance(a2), Catch::Matchers::WithinAbs(std::sqrt(2), 1e-6));
    CHECK_THAT(a1.distance(a3), Catch::Matchers::WithinAbs(std::sqrt(2), 1e-6));
    CHECK_THAT(a1.distance(a4), Catch::Matchers::WithinAbs(std::sqrt(2), 1e-6));
    CHECK(a1.distance(a5) == 0);
}

TEST_CASE("operators") {
//*** ATOMS ***//
    Atom a1 = Atom({3, 0, 5}, 2, "He", "", 3);
    Atom a2 = a1;
    REQUIRE(a1 == a2);
    REQUIRE(a1.equals_content(a2));

    a2 = Atom({0, 4, 1}, 2, "He", "", 2);
    REQUIRE(a1 != a2);
    REQUIRE(!a1.equals_content(a2));
    REQUIRE(a2 < a1);

//*** HETATOMS ***//
    Water w1 = Water({3, 0, 5}, 2, "He", "", 3);
    Water w2 = w1;
    REQUIRE(w1 == w2);

    w2 = Atom({0, 4, 1}, 2, "He", "", 2);
    REQUIRE(w1 != w2);
    REQUIRE(w2 < w1);
}