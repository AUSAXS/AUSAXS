#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <data/Molecule.h>
#include <data/Body.h>
#include <data/record/Water.h>
#include <data/record/Record.h>
#include <utility/Console.h>
#include <settings/All.h>

#include <vector>
#include <string>
#include <fstream>
#include <filesystem>
#include <iostream>

using namespace ausaxs;
using namespace data;
using namespace data::record;

TEST_CASE("PDBReader::read") {
    settings::molecule::center = false;
    settings::general::verbose = false;

    io::File path("temp/io/temp.pdb");
    path.create();

    std::ofstream pdb_file(path);
    pdb_file << "ATOM      1  CB  ARG A 129         2.1     3.2     4.3  0.50 42.04           C " << std::endl;
    pdb_file << "ATOM      2  CB  ARG A 129         3.2     4.3     5.4  0.50 42.04           C " << std::endl;
    pdb_file << "TER       3      ARG A 129                                                     " << std::endl;
    pdb_file << "HETATM    4  O   HOH A 130      30.117  29.049  34.879  0.94 34.19           O " << std::endl;
    pdb_file << "HETATM    5  O   HOH A 131      31.117  30.049  35.879  0.94 34.19           O " << std::endl;
    pdb_file.close();

    // check PDB io
    std::unique_ptr<Molecule> protein = std::make_unique<Molecule>(path);
    path = "temp/io/temp2.pdb";
    protein->save(path);
    protein = std::make_unique<Molecule>(path);
    std::vector<Atom> atoms = protein->get_atoms();
    Atom a = atoms[0];

    if (atoms.size() == 0) {
        REQUIRE(false);
        return;
    }

    // the idea is that we have now loaded the hardcoded strings above, saved them, and loaded them again. 
    // we now compare the loaded values with the expected.
    CHECK(a.serial == 1);
    CHECK(a.coords.x() == 2.1);
    CHECK(a.coords.y() == 3.2);
    CHECK(a.coords.z() == 4.3);
    CHECK(a.occupancy == 0.50);
    CHECK(a.element == constants::atom_t::C);
    CHECK(a.resName == "ARG");
}

TEST_CASE("PDBWriter: writing multifile pdb") {
    settings::molecule::implicit_hydrogens = false;

    std::vector<Atom> atoms(101000);
    for (unsigned int i = 0; i < atoms.size(); i++) {
        atoms[i] = Atom({1,2,3}, 1, constants::atom_t::C, "LYS", i);
    }
    std::vector<Water> waters(100);
    for (unsigned int i = 0; i < waters.size(); i++) {
        waters[i] = Water({1,2,3}, 1, constants::atom_t::O, "HOH", i);
    }

    Molecule protein(atoms, waters);
    protein.save("temp/io/temp_multifile.pdb");
    
    REQUIRE(std::filesystem::exists("temp/io/temp_multifile_part1.pdb"));
    REQUIRE(std::filesystem::exists("temp/io/temp_multifile_part2.pdb"));

    // first file
    Molecule protein2("temp/io/temp_multifile_part1.pdb");
    REQUIRE(protein2.get_body(0).size_atom() == 100000);
    REQUIRE(protein2.get_body(0).get_atoms().back().serial == 99999);

    Molecule protein3("temp/io/temp_multifile_part2.pdb");
    REQUIRE(protein3.get_body(0).size_atom() == 1000);
    for (int i = 0; i < 1000; i++) {
        REQUIRE(protein3.get_body(0).get_atom(i).serial == i);
    }
    REQUIRE(protein3.get_body(0).get_file().terminate.serial == 1000);
    REQUIRE(protein3.size_water() == 100);
    for (int i = 0; i < 100; i++) {
        REQUIRE(protein3.get_water(i).serial == i+1001);
    }
}

TEST_CASE("PDBReader: can_parse_hydrogens") {
    std::vector<std::string> val = {"ATOM      1  N   VAL     1      -3.299   8.066 -11.443  1.00  0.00           N",
                                    "ATOM      2  H   VAL     1      -3.411   8.677 -12.239  1.00  0.00           H",
                                    "ATOM      3  CA  VAL     1      -3.085   6.673 -11.780  1.00  0.00           C",
                                    "ATOM      4  HA  VAL     1      -3.328   6.080 -10.899  1.00  0.00           H",
                                    "ATOM      5  CB  VAL     1      -3.927   6.165 -12.947  1.00  0.00           C",
                                    "ATOM      6  HB  VAL     1      -3.774   6.930 -13.708  1.00  0.00           H",
                                    "ATOM      7  CG1 VAL     1      -3.577   4.780 -13.486  1.00  0.00           C",
                                    "ATOM      8 HG11 VAL     1      -3.508   4.011 -12.716  1.00  0.00           H",
                                    "ATOM      9 HG12 VAL     1      -4.289   4.438 -14.237  1.00  0.00           H",
                                    "ATOM     10 HG13 VAL     1      -2.612   4.760 -13.992  1.00  0.00           H",
                                    "ATOM     11  CG2 VAL     1      -5.370   6.065 -12.463  1.00  0.00           C",
                                    "ATOM     12 HG21 VAL     1      -5.876   7.021 -12.328  1.00  0.00           H",
                                    "ATOM     13 HG22 VAL     1      -6.043   5.567 -13.162  1.00  0.00           H",
                                    "ATOM     14 HG23 VAL     1      -5.354   5.510 -11.524  1.00  0.00           H",
                                    "ATOM     15  C   VAL     1      -1.621   6.436 -12.123  1.00  0.00           C",
                                    "ATOM     16  O   VAL     1      -1.200   7.045 -13.104  1.00  0.00           O",
    };

    Atom atom;
    atom.parse_pdb(val[7]);
    REQUIRE(atom.name == "HG11");
    REQUIRE(atom.element == constants::atom_t::H);
    REQUIRE_THAT(atom.coords.x(), Catch::Matchers::WithinAbs(-3.508, 1e-6));
    REQUIRE_THAT(atom.coords.y(), Catch::Matchers::WithinAbs(4.011, 1e-6));
    REQUIRE_THAT(atom.coords.z(), Catch::Matchers::WithinAbs(-12.716, 1e-6));

    atom.parse_pdb(val[8]);
    REQUIRE(atom.name == "HG12");
    REQUIRE(atom.element == constants::atom_t::H);
    REQUIRE_THAT(atom.coords.x(), Catch::Matchers::WithinAbs(-4.289, 1e-6));
    REQUIRE_THAT(atom.coords.y(), Catch::Matchers::WithinAbs(4.438, 1e-6));
    REQUIRE_THAT(atom.coords.z(), Catch::Matchers::WithinAbs(-14.237, 1e-6));

    atom.parse_pdb(val[9]);
    REQUIRE(atom.name == "HG13");
    REQUIRE(atom.element == constants::atom_t::H);
    REQUIRE_THAT(atom.coords.x(), Catch::Matchers::WithinAbs(-2.612, 1e-6));
    REQUIRE_THAT(atom.coords.y(), Catch::Matchers::WithinAbs(4.760, 1e-6));
    REQUIRE_THAT(atom.coords.z(), Catch::Matchers::WithinAbs(-13.992, 1e-6));
}