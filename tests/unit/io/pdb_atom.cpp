#include <catch2/catch_test_macros.hpp>

#include <io/pdb/PDBAtom.h>
#include <io/pdb/PDBWater.h>
#include <constants/Constants.h>

using namespace ausaxs;
using namespace ausaxs::io::pdb;

TEST_CASE("PDBAtom::PDBAtom") {
    SECTION("Vector3<double>, double, std::string&, std::string&, int") {
        PDBAtom a1({3, 0, 5}, 2, constants::atom_t::He, "LYS", 3);

        CHECK(a1.serial == 3);
        CHECK(a1.resName == "LYS");
        CHECK(a1.coordinates() == Vector3<double>{3, 0, 5});
        CHECK(a1.occupancy == 2);
        CHECK(a1.element == constants::atom_t::He);
        CHECK(a1.is_water() == false);
    }

    SECTION("int, std::string&, std::string&, std::string&, int, std::string&, Vector3<double>, double, double, std::string&, std::string&") {
        PDBAtom a1(15, "CA", "altLoc", "GLY", 'X', 3, "iCode", Vector3<double>{0, 1, 2}, 2.5, 3.5, constants::atom_t::He, "2-");

        CHECK(a1.serial == 15);
        CHECK(a1.name == "CA");
        CHECK(a1.altLoc == "altLoc");
        CHECK(a1.resName == "GLY");
        CHECK(a1.chainID == 'X');
        CHECK(a1.resSeq == 3);
        CHECK(a1.iCode == "iCode");
        CHECK(a1.coordinates() == Vector3<double>{0, 1, 2});
        CHECK(a1.occupancy == 2.5);
        CHECK(a1.tempFactor == 3.5);
        CHECK(a1.element == constants::atom_t::He);
        CHECK(a1.charge == "2-");
        CHECK(a1.is_water() == false);
    }    
}

TEST_CASE("PDBAtom::get_type") {
    SECTION("PDBAtom") {
        PDBAtom a1({3, 0, 5}, 2, constants::atom_t::He, "resName", 3);
        CHECK(a1.get_type() == RecordType::ATOM);
    }
    SECTION("PDBWater") {
        PDBWater w1 = PDBWater::create_new_water(Vector3<double>{1, 2, 3});
        CHECK(w1.get_type() == RecordType::WATER);
    }
}

TEST_CASE("PDBAtom::parse_pdb") {
    SECTION("atom") {
        std::string line = "ATOM      1  N   GLY A   1       1.000   2.000   3.000  1.00  0.00           N  ";
        PDBAtom a1; a1.parse_pdb(line);
        CHECK(a1.serial == 1);
        CHECK(a1.name == "N");
        CHECK(a1.altLoc == " ");
        CHECK(a1.resName == "GLY");
        CHECK(a1.chainID == 'A');
        CHECK(a1.resSeq == 1);
        CHECK(a1.iCode == " ");
        CHECK(a1.coordinates() == Vector3<double>{1, 2, 3});
        CHECK(a1.occupancy == 1);
        CHECK(a1.tempFactor == 0);
        CHECK(a1.element == constants::atom_t::N);
        CHECK(a1.charge == "  ");
        CHECK(a1.is_water() == false);
    }

    SECTION("hetatm") {
        std::string line = "HETATM    2  C   LYS B   2       5.000   4.000   2.000  0.50  0.50           C  ";
        PDBAtom a1; a1.parse_pdb(line);
        CHECK(a1.serial == 2);
        CHECK(a1.name == "C");
        CHECK(a1.altLoc == " ");
        CHECK(a1.resName == "LYS");
        CHECK(a1.chainID == 'B');
        CHECK(a1.resSeq == 2);
        CHECK(a1.iCode == " ");
        CHECK(a1.coordinates() == Vector3<double>{5, 4, 2});
        CHECK(a1.occupancy == 0.5);
        CHECK(a1.tempFactor == 0.5);
        CHECK(a1.element == constants::atom_t::C);
        CHECK(a1.charge == "  ");
        CHECK(a1.is_water() == false);
    }

    SECTION("custom") {
        std::string line = "ATOM      1  CB  ARG A 129         2.1     3.2     4.3  0.50 42.04           C ";
        PDBAtom a; a.parse_pdb(line);
        CHECK(a.serial == 1);
        CHECK(a.name == "CB");
        CHECK(a.altLoc == " ");
        CHECK(a.resName == "ARG");
        CHECK(a.chainID == 'A');
        CHECK(a.resSeq == 129);
        CHECK(a.iCode == " ");
        CHECK(a.coordinates() == Vector3<double>{2.1, 3.2, 4.3});
        CHECK(a.occupancy == 0.5);
        CHECK(a.tempFactor == 42.04);
        CHECK(a.element == constants::atom_t::C);
        CHECK(a.charge == "  ");
        CHECK(a.get_recName() == "ATOM  ");
    }
}

TEST_CASE("PDBAtom::as_pdb") {
    SECTION("parse and as_pdb") {
        PDBAtom a1; a1.parse_pdb("ATOM      1  N   GLY A   1     1.00000 2.00000 3.000001.00000.0000           N  ");
        std::string res = a1.as_pdb();
        PDBAtom a2; a2.parse_pdb(res);
        CHECK(a1.equals_content(a2));
    }
}

TEST_CASE("PDBAtom::coordinates") {
    PDBAtom a1({1, 2, 3}, 1, constants::atom_t::N, "GLY", 1);
    CHECK(a1.coordinates() == Vector3<double>{1, 2, 3});
    
    a1.coordinates() = Vector3<double>{4, 5, 6};
    CHECK(a1.coordinates() == Vector3<double>{4, 5, 6});
}

TEST_CASE("PDBAtom::get_recName") {
    PDBAtom a1;
    CHECK(a1.get_recName() == "ATOM  ");
}

TEST_CASE("PDBAtom::get_mass") {
    SECTION("H") {
        PDBAtom a1;
        a1.set_element("H");
        a1.resName = "GLY";
        a1.name = "H";
        CHECK(a1.get_mass() == constants::mass::get_mass(constants::atom_t::H));
    }

    SECTION("C") {
        PDBAtom a1;
        a1.set_element("C");
        a1.resName = "GLY";
        a1.name = "CA"; // CA has 2 H attached
        CHECK(a1.get_mass() == constants::mass::get_mass(constants::atom_t::C) + 2*constants::mass::get_mass(constants::atom_t::H));
    }
}

TEST_CASE("PDBAtom::Z") {
    SECTION("H") {
        PDBAtom a1;
        a1.set_element("H");
        CHECK(a1.Z() == 1);
    }

    SECTION("C") {
        PDBAtom a1;
        a1.set_element("C");
        CHECK(a1.Z() == 6);
    }
}

TEST_CASE("PDBAtom::operator<") {
    PDBAtom a1;
    PDBAtom a2;
    a1.serial = 1;
    a2.serial = 2;
    CHECK(a1 < a2);
    
    a2.serial = 1;
    CHECK(!(a1 < a2));
}

TEST_CASE("PDBAtom::operator==") {
    PDBAtom a1;
    PDBAtom a2;
    CHECK(!(a1 == a2));

    PDBAtom a3;
    CHECK(!(a1 == a3));
}

TEST_CASE("PDBAtom::equals_content") {
    SECTION("simple") {
        PDBAtom a1;
        PDBAtom a2;
        CHECK(a1.equals_content(a2));

        a1.serial = 5;
        CHECK(!a1.equals_content(a2));

        a2.serial = 5;
        CHECK(a1.equals_content(a2));
    }

    SECTION("complex") {
        std::string line = "ATOM      1  N   GLY A   1      11.000  12.000  13.000  1.00  0.00           N  ";
        PDBAtom a1; a1.parse_pdb(line);
        PDBAtom a2; a2.parse_pdb(line);
        CHECK(a1.equals_content(a2));

        a1.serial = 5;
        CHECK(!a1.equals_content(a2));

        a2.serial = 5;
        CHECK(a1.equals_content(a2));
    }
}

TEST_CASE("PDBAtom::set_element") {
    PDBAtom a1;
    a1.set_element("C");
    CHECK(a1.element == constants::atom_t::C);
    
    a1.set_element(constants::atom_t::N);
    CHECK(a1.element == constants::atom_t::N);
}
