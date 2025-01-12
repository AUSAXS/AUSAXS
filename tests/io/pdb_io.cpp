#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <data/Molecule.h>
#include <data/Body.h>
#include <io/pdb/PDBStructure.h>
#include <io/detail/PDBReader.h>
#include <io/detail/PDBWriter.h>
#include <utility/Console.h>
#include <settings/All.h>

#include <vector>
#include <string>
#include <fstream>
#include <iostream>

using namespace ausaxs;
using namespace ausaxs::data;
using namespace ausaxs::io::pdb;

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
    auto protein = io::detail::pdb::read(path);
    path = "temp/io/temp2.pdb";
    io::detail::pdb::write(protein, path);
    protein = io::detail::pdb::read(path);

    // the idea is that we have now loaded the hardcoded strings above, saved them, and loaded them again. 
    // we now compare the loaded values with the expected.
    REQUIRE(protein.atoms.size() == 2);
    CHECK(protein.atoms[0].serial == 1);
    CHECK(protein.atoms[0].coords.x() == 2.1);
    CHECK(protein.atoms[0].coords.y() == 3.2);
    CHECK(protein.atoms[0].coords.z() == 4.3);
    CHECK(protein.atoms[0].occupancy == 0.50);
    CHECK(protein.atoms[0].element == constants::atom_t::C);
    CHECK(protein.atoms[0].resName == "ARG");

    REQUIRE(protein.waters.size() == 2);
    CHECK(protein.waters[0].serial == 4);
    CHECK(protein.waters[0].coords.x() == 30.117);
    CHECK(protein.waters[0].coords.y() == 29.049);
    CHECK(protein.waters[0].coords.z() == 34.879);
    CHECK(protein.waters[0].occupancy == 0.94);
    CHECK(protein.waters[0].element == constants::atom_t::O);
    CHECK(protein.waters[0].resName == "HOH");
}

TEST_CASE("PDBReader: add_implicit_hydrogens", "[files]") {
    settings::molecule::implicit_hydrogens = true;
    std::vector<PDBAtom> atoms = {
        PDBAtom(1, "N",  "", "LYS", 'A', 1, "", Vector3<double>(0, 0, 0), 1, 0, constants::atom_t::N, "0"),
        PDBAtom(2, "CA", "", "LYS", 'A', 1, "", Vector3<double>(0, 0, 0), 1, 0, constants::atom_t::C, "0"),
        PDBAtom(3, "C",  "", "LYS", 'A', 1, "", Vector3<double>(0, 0, 0), 1, 0, constants::atom_t::C, "0"),
        PDBAtom(4, "O",  "", "LYS", 'A', 1, "", Vector3<double>(0, 0, 0), 1, 0, constants::atom_t::O, "0"),
        PDBAtom(5, "CB", "", "LYS", 'A', 1, "", Vector3<double>(0, 0, 0), 1, 0, constants::atom_t::C, "0"),
        PDBAtom(6, "CG", "", "LYS", 'A', 1, "", Vector3<double>(0, 0, 0), 1, 0, constants::atom_t::C, "0"),
        PDBAtom(7, "CD", "", "LYS", 'A', 1, "", Vector3<double>(0, 0, 0), 1, 0, constants::atom_t::C, "0"),
        PDBAtom(8, "CE", "", "LYS", 'A', 1, "", Vector3<double>(0, 0, 0), 1, 0, constants::atom_t::C, "0"),
        PDBAtom(9, "NZ", "", "LYS", 'A', 1, "", Vector3<double>(0, 0, 0), 1, 0, constants::atom_t::N, "0"),
    };
    io::pdb::PDBStructure protein({atoms, {}});

    for (auto a : protein.get_atoms()) {
        CHECK(a.effective_charge == constants::charge::nuclear::get_charge(a.get_element()));
        CHECK(a.get_atomic_group() == constants::atomic_group_t::unknown);
    }

    protein.add_implicit_hydrogens();
    atoms = protein.get_atoms();
    CHECK(atoms[0].effective_charge == constants::charge::nuclear::get_charge(atoms[0].get_element()) + 1);
    CHECK(atoms[0].get_atomic_group() == constants::atomic_group_t::NH);

    CHECK(atoms[1].effective_charge == constants::charge::nuclear::get_charge(atoms[1].get_element()) + 1);
    CHECK(atoms[1].get_atomic_group() == constants::atomic_group_t::CH);

    CHECK(atoms[2].effective_charge == constants::charge::nuclear::get_charge(atoms[2].get_element()) + 0);
    CHECK(atoms[2].get_atomic_group() == constants::atomic_group_t::unknown);

    CHECK(atoms[3].effective_charge == constants::charge::nuclear::get_charge(atoms[3].get_element()) + 0);
    CHECK(atoms[3].get_atomic_group() == constants::atomic_group_t::unknown);

    CHECK(atoms[4].effective_charge == constants::charge::nuclear::get_charge(atoms[4].get_element()) + 2);
    CHECK(atoms[4].get_atomic_group() == constants::atomic_group_t::CH2);

    CHECK(atoms[5].effective_charge == constants::charge::nuclear::get_charge(atoms[5].get_element()) + 2);
    CHECK(atoms[5].get_atomic_group() == constants::atomic_group_t::CH2);

    CHECK(atoms[6].effective_charge == constants::charge::nuclear::get_charge(atoms[6].get_element()) + 2);
    CHECK(atoms[6].get_atomic_group() == constants::atomic_group_t::CH2);

    CHECK(atoms[7].effective_charge == constants::charge::nuclear::get_charge(atoms[7].get_element()) + 2);
    CHECK(atoms[7].get_atomic_group() == constants::atomic_group_t::CH2);

    CHECK(atoms[8].effective_charge == constants::charge::nuclear::get_charge(atoms[8].get_element()) + 3);
    CHECK(atoms[8].get_atomic_group() == constants::atomic_group_t::NH3);
}

TEST_CASE("PDBWriter: writing multifile pdb") {
    settings::molecule::implicit_hydrogens = false;

    std::vector<AtomFF> atoms(101000);
    for (int i = 0; i < static_cast<int>(atoms.size()); ++i) {
        atoms[i] = AtomFF({1, 2, 3}, form_factor::form_factor_t::C);
    }
    std::vector<Water> waters(100);
    for (int i = 0; i < static_cast<int>(waters.size()); ++i) {
        waters[i] = Water({1, 2, 3});
    }

    Molecule protein({{atoms, waters}});
    protein.save("temp/io/temp_multifile.pdb");
    
    REQUIRE(io::File("temp/io/temp_multifile_part1.pdb").exists());
    REQUIRE(io::File("temp/io/temp_multifile_part2.pdb").exists());

    // first file
    Molecule protein2("temp/io/temp_multifile_part1.pdb");
    REQUIRE(protein2.get_body(0).size_atom() == 100000);

    Molecule protein3("temp/io/temp_multifile_part2.pdb");
    REQUIRE(protein3.get_body(0).size_atom() == 1000);
    REQUIRE(protein3.size_water() == 100);

    Molecule protein4("temp/io/temp_multifile.pdb");
    REQUIRE(protein.size_body() == protein4.size_body());
    for (unsigned int i = 0; i < protein.get_bodies().size(); ++i) {
        REQUIRE(protein.get_body(i).equals_content(protein4.get_body(i)));
    }
}

TEST_CASE("PDBReader: can_parse_hydrogens") {
    std::vector<std::string> val = {
        "ATOM      1  N   VAL     1      -3.299   8.066 -11.443  1.00  0.00           N",
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

    PDBAtom atom;
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