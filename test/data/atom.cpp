#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>

#include <data/Atom.h>
#include <data/Water.h>
#include <settings/All.h>
#include <utility/Constants.h>

TEST_CASE("Atom::Atom") {
    SECTION("Vector3<double>, double, std::string&, std::string&, int") {
        Atom a1 = Atom({3, 0, 5}, 2, "He", "LYS", 3);

        CHECK(a1.get_serial() == 3);
        CHECK(a1.get_resName() == "LYS");
        CHECK(a1.get_coordinates() == Vector3<double>({3, 0, 5}));
        CHECK(a1.get_occupancy() == 2);
        CHECK(a1.get_element() == "He");
        CHECK(a1.is_water() == false);
    }

    SECTION("int, std::string&, std::string&, std::string&, int, std::string&, Vector3<double>, double, double, std::string&, std::string&") {
        Atom a1 = Atom(15, "CA", "altLoc", "GLY", "chainID", 3, "iCode", Vector3<double>({0, 1, 2}), 2.5, 3.5, "He", "2-");

        CHECK(a1.get_serial() == 15);
        CHECK(a1.get_name() == "CA");
        CHECK(a1.get_altLoc() == "altLoc");
        CHECK(a1.get_resName() == "GLY");
        CHECK(a1.get_chainID() == "chainID");
        CHECK(a1.get_resSeq() == 3);
        CHECK(a1.get_iCode() == "iCode");
        CHECK(a1.get_coordinates() == Vector3({0, 1, 2}));
        CHECK(a1.get_occupancy() == 2.5);
        CHECK(a1.get_tempFactor() == 3.5);
        CHECK(a1.get_element() == "He");
        CHECK(a1.get_charge() == "2-");
        CHECK(a1.is_water() == false);
    }    
}

TEST_CASE("Atom::get_type") {
    SECTION("Atom") {
        Atom a1 = Atom({3, 0, 5}, 2, "He", "resName", 3);
        CHECK(a1.get_type() == RecordType::ATOM);
    }
    SECTION("Water") {
        Water w1 = Water::create_new_water(Vector3<double>({1, 2, 3}));
        CHECK(w1.get_type() == RecordType::WATER);
    }
}

TEST_CASE("Atom::parse_pdb") {
    SECTION("atom") {
        std::string line = "ATOM      1  N   GLY A   1       1.000   2.000   3.000  1.00  0.00           N  ";
        Atom a1; a1.parse_pdb(line);
        CHECK(a1.get_serial() == 1);
        CHECK(a1.get_name() == "N");
        CHECK(a1.get_altLoc() == " ");
        CHECK(a1.get_resName() == "GLY");
        CHECK(a1.get_chainID() == "A");
        CHECK(a1.get_resSeq() == 1);
        CHECK(a1.get_iCode() == " ");
        CHECK(a1.get_coordinates() == Vector3({1, 2, 3}));
        CHECK(a1.get_occupancy() == 1);
        CHECK(a1.get_tempFactor() == 0);
        CHECK(a1.get_element() == "N");
        CHECK(a1.get_charge() == "  ");
        CHECK(a1.is_water() == false);
    }

    SECTION("hetatm") {
        std::string line = "HETATM    2  C   LYS B   2       5.000   4.000   2.000  0.50  0.50           C  ";
        Atom a1; a1.parse_pdb(line);
        CHECK(a1.get_serial() == 2);
        CHECK(a1.get_name() == "C");
        CHECK(a1.get_altLoc() == " ");
        CHECK(a1.get_resName() == "LYS");
        CHECK(a1.get_chainID() == "B");
        CHECK(a1.get_resSeq() == 2);
        CHECK(a1.get_iCode() == " ");
        CHECK(a1.get_coordinates() == Vector3({5, 4, 2}));
        CHECK(a1.get_occupancy() == 0.5);
        CHECK(a1.get_tempFactor() == 0.5);
        CHECK(a1.get_element() == "C");
        CHECK(a1.get_charge() == "  ");
        CHECK(a1.is_water() == false);
    }

    SECTION("custom") {
        std::string line = "ATOM      1  CB  ARG A 129         2.1     3.2     4.3  0.50 42.04           C ";
        Atom a; a.parse_pdb(line);
        CHECK(a.get_serial() == 1);
        CHECK(a.get_name() == "CB");
        CHECK(a.get_altLoc() == " "); // spaces are only removed from number strings
        CHECK(a.get_resName() == "ARG");
        CHECK(a.get_chainID() == "A");
        CHECK(a.get_resSeq() == 129);
        CHECK(a.get_iCode() == " "); // same
        CHECK(a.get_coordinates() == Vector3({2.1, 3.2, 4.3}));
        CHECK(a.get_occupancy() == 0.5);
        CHECK(a.get_tempFactor() == 42.04);
        CHECK(a.get_element() == "C");
        CHECK(a.get_charge() == "  "); // same
        CHECK(a.get_recName() == "ATOM  ");
    }
}

TEST_CASE("Atom::as_pdb") {
    SECTION("atom") {
        Atom a1 = Atom({1, 2, 3}, 1, "N", "GLY", 1);
        std::string line = a1.as_pdb().substr(0, 80);
        std::string res = "ATOM      1      GLY    -1     1.00000 2.00000 3.000001.0000-1.000           N  ";
        CHECK(line == res);
    }

    SECTION("parse and as_pdb") {
        std::string line = "ATOM      1  N   GLY A   1     1.00000 2.00000 3.000001.00000.0000           N  ";
        Atom a1; a1.parse_pdb(line);
        std::string res = a1.as_pdb().substr(0, 80);
        CHECK(res == line);
    }
}

TEST_CASE("Atom::distance") {
    SECTION("simple") {
        Atom a1 = Atom({0, 0, 0}, 1, "N", "GLY", 1);
        Atom a2 = Atom({1, 0, 0}, 2, "C", "GLY", 1);
        CHECK(a1.distance(a2) == 1);
    }

    SECTION("complex") {
        Atom a1 = Atom({0, 0, 0}, 1, "N", "GLY", 1);
        Atom a2 = Atom({1, 0, 0}, 2, "C", "GLY", 1);
        Atom a3 = Atom({1, 1, 0}, 3, "O", "GLY", 1);
        CHECK(a1.distance(a2) == 1);
        CHECK(a1.distance(a3) == std::sqrt(2));
    }
}

TEST_CASE("Atom::translate") {
    SECTION("simple") {
        Atom a1 = Atom({0, 0, 0}, 1, "N", "GLY", 1);
        a1.translate({1, 2, 3});
        CHECK(a1.get_coordinates() == Vector3({1, 2, 3}));
    }

    SECTION("complex") {
        Atom a1 = Atom({0, 0, 0}, 1, "N", "GLY", 1);
        a1.translate({1, 2, 3});
        CHECK(a1.get_coordinates() == Vector3({1, 2, 3}));
        a1.translate({1, 2, 3});
        CHECK(a1.get_coordinates() == Vector3({2, 4, 6}));
    }
}

TEST_CASE("Atom::set_coordinates") {
    Atom a1;
    a1.set_coordinates({1, 2, 3});
    CHECK(a1.get_coordinates() == Vector3({1, 2, 3}));
}

TEST_CASE("Atom::set_x") {
    Atom a1;
    a1.set_x(1);
    CHECK(a1.get_coordinates() == Vector3({1, 0, 0}));
}

TEST_CASE("Atom::set_y") {
    Atom a1;
    a1.set_y(1);
    CHECK(a1.get_coordinates() == Vector3({0, 1, 0}));
}

TEST_CASE("Atom::set_z") {
    Atom a1;
    a1.set_z(1);
    CHECK(a1.get_coordinates() == Vector3({0, 0, 1}));
}

TEST_CASE("Atom::set_occupancy") {
    Atom a1;
    a1.set_occupancy(2);
    CHECK(a1.get_occupancy() == 2);
}

TEST_CASE("Atom::set_tempFactor") {
    Atom a1;
    a1.set_tempFactor(2);
    CHECK(a1.get_tempFactor() == 2);
}

TEST_CASE("Atom::set_altLoc") {
    Atom a1;
    a1.set_altLoc("Z");
    CHECK(a1.get_altLoc() == "Z");
}

TEST_CASE("Atom::set_serial") {
    Atom a1;
    a1.set_serial(2);
    CHECK(a1.get_serial() == 2);
}

TEST_CASE("Atom::set_resSeq") {
    Atom a1;
    a1.set_resSeq(2);
    CHECK(a1.get_resSeq() == 2);
}

TEST_CASE("Atom::set_iCode") {
    Atom a1;
    a1.set_iCode("Z");
    CHECK(a1.get_iCode() == "Z");
}

TEST_CASE("Atom::set_chainID") {
    Atom a1;
    a1.set_chainID("Z");
    CHECK(a1.get_chainID() == "Z");
}

TEST_CASE("Atom::set_element") {
    Atom a1;
    a1.set_element("He");
    CHECK(a1.get_element() == "He");
}

TEST_CASE("Atom::set_charge") {
    Atom a1;
    a1.set_charge("2-");
    CHECK(a1.get_charge() == "2-");
}

TEST_CASE("Atom::set_resName") {
    Atom a1;
    a1.set_resName("resName");
    CHECK(a1.get_resName() == "resName");
}

TEST_CASE("Atom::set_name") {
    Atom a1;
    a1.set_name("name");
    CHECK(a1.get_name() == "name");
}

TEST_CASE("Atom::set_effective_charge") {
    Atom a1;
    a1.set_effective_charge(2);
    CHECK(a1.get_effective_charge() == 2);
}

TEST_CASE("Atom::get_recName") {
    Atom a1;
    CHECK(a1.get_recName() == "ATOM  ");
}

TEST_CASE("Atom::get_mass") {
    SECTION("H") {
        Atom a1;
        a1.set_element("H");
        a1.set_resName("GLY");
        a1.set_name("H");
        CHECK(a1.get_mass() == constants::mass::atomic.get("H"));
    }

    SECTION("C") {
        Atom a1;
        a1.set_element("C");
        a1.set_resName("GLY");
        a1.set_name("CA"); // CA has 2 H attached
        CHECK(a1.get_mass() == constants::mass::atomic.get("C") + 2*constants::mass::atomic.get("H"));
    }
}

TEST_CASE("Atom::Z") {
    SECTION("H") {
        Atom a1;
        a1.set_element("H");
        CHECK(a1.Z() == 1);
        CHECK(a1.get_absolute_charge() == a1.Z());
    }

    SECTION("C") {
        Atom a1;
        a1.set_element("C");
        CHECK(a1.Z() == 6);
        CHECK(a1.get_absolute_charge() == a1.Z());
    }
}

// compares by serial
TEST_CASE("Atom::operator<") {
    Atom a1;
    Atom a2;
    a1.set_serial(1);
    a2.set_serial(2);
    CHECK(a1 < a2);
    
    a2.set_serial(1);
    CHECK(!(a1 < a2));
}

// compares by uid
TEST_CASE("Atom::operator==") {
    Atom a1;
    Atom a2;
    CHECK(!(a1 == a2));

    Atom a3;
    CHECK(!(a1 == a3));
}

TEST_CASE("Atom::equals_content") {
    SECTION("simple") {
        Atom a1;
        Atom a2;
        CHECK(a1.equals_content(a2));

        a1.set_serial(5);
        CHECK(!a1.equals_content(a2));

        a2.set_serial(5);
        CHECK(a1.equals_content(a2));
    }

    SECTION("complex") {
        std::string line = "ATOM      1  N   GLY A   1      11.000  12.000  13.000  1.00  0.00           N  ";
        Atom a1; a1.parse_pdb(line);
        Atom a2; a2.parse_pdb(line);
        CHECK(a1.equals_content(a2));

        a1.set_serial(5);
        CHECK(!a1.equals_content(a2));

        a2.set_serial(5);
        CHECK(a1.equals_content(a2));
    }
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