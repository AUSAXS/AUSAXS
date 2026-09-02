#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <data/Molecule.h>
#include <data/Body.h>
#include <io/pdb/PDBStructure.h>
#include <io/detail/structure/PDBReader.h>
#include <io/detail/structure/PDBWriter.h>
#include <utility/Console.h>
#include <settings/All.h>

#include <support/temp_file.h>

#include <vector>
#include <string>
#include <iostream>

using namespace ausaxs;
using namespace ausaxs::data;
using namespace ausaxs::io::pdb;

TEST_CASE("PDBReader::read") {
    settings::molecule::center = false;
    settings::general::verbose = false;

    std::string content =
        "ATOM      1  CB  ARG A 129         2.1     3.2     4.3  0.50 42.04           C \n"
        "ATOM      2  CB  ARG A 129         3.2     4.3     5.4  0.50 42.04           C \n"
        "TER       3      ARG A 129                                                     \n"
        "HETATM    4  O   HOH A 130      30.117  29.049  34.879  0.94 34.19           O \n"
        "HETATM    5  O   HOH A 131      31.117  30.049  35.879  0.94 34.19           O \n";
    test::TempFile path(".pdb", content);

    // check PDB io
    auto protein = io::detail::pdb::read(path);
    test::TempFile path2(".pdb");
    io::detail::pdb::write(protein, path2);
    protein = io::detail::pdb::read(path2);

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

TEST_CASE("PDBReader: add_implicit_hydrogens") {
    auto generate_molecule = [] () {
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
        return io::pdb::PDBStructure({atoms, {}});
    };

    SECTION("enabled") {
        auto protein = generate_molecule();
        protein.add_implicit_hydrogens();
        auto& atoms = protein.atoms;

        CHECK(atoms[0].effective_charge == constants::charge::get_ff_charge(atoms[0].get_form_factor_type()) + 1);
        CHECK(atoms[0].atomic_group == constants::atomic_group_t::NH);

        CHECK(atoms[1].effective_charge == constants::charge::get_ff_charge(atoms[1].get_form_factor_type()) + 1);
        CHECK(atoms[1].atomic_group == constants::atomic_group_t::CH);

        CHECK(atoms[2].effective_charge == constants::charge::get_ff_charge(atoms[2].get_form_factor_type()) + 0);
        CHECK(atoms[2].atomic_group == constants::atomic_group_t::unknown);

        CHECK(atoms[3].effective_charge == constants::charge::get_ff_charge(atoms[3].get_form_factor_type()) + 0);
        CHECK(atoms[3].atomic_group == constants::atomic_group_t::unknown);

        CHECK(atoms[4].effective_charge == constants::charge::get_ff_charge(atoms[4].get_form_factor_type()) + 2);
        CHECK(atoms[4].atomic_group == constants::atomic_group_t::CH2);

        CHECK(atoms[5].effective_charge == constants::charge::get_ff_charge(atoms[5].get_form_factor_type()) + 2);
        CHECK(atoms[5].atomic_group == constants::atomic_group_t::CH2);

        CHECK(atoms[6].effective_charge == constants::charge::get_ff_charge(atoms[6].get_form_factor_type()) + 2);
        CHECK(atoms[6].atomic_group == constants::atomic_group_t::CH2);

        CHECK(atoms[7].effective_charge == constants::charge::get_ff_charge(atoms[7].get_form_factor_type()) + 2);
        CHECK(atoms[7].atomic_group == constants::atomic_group_t::CH2);

        CHECK(atoms[8].effective_charge == constants::charge::get_ff_charge(atoms[8].get_form_factor_type()) + 3);
        CHECK(atoms[8].atomic_group == constants::atomic_group_t::NH3);
    }

    SECTION("disabled") {
        auto protein = generate_molecule();

        for (auto a : protein.atoms) {
            CHECK(a.effective_charge == constants::charge::get_ff_charge(a.get_form_factor_type()));
            CHECK(a.atomic_group == constants::atomic_group_t::unknown);
        }
    }
}

TEST_CASE("PDBWriter: multi-chain structures keep their chains apart") {
    settings::general::verbose = false;
    settings::molecule::center = false;
    settings::molecule::implicit_hydrogens = false;
    settings::molecule::store_calpha = true;
    settings::molecule::store_residue_seq = true;
    settings::molecule::store_chain_id = true;

    // Molecule loads an entire file into a single Body, so the chain boundaries only survive as metadata. The residue ids restart at every chain, so writing
    // the body under one chain identifier would make the residues of the second chain indistinguishable from those of the first.
    std::string content =
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N \n"
        "ATOM      2  CA  ALA A   1       1.000   0.000   0.000  1.00  0.00           C \n"
        "ATOM      3  C   ALA A   1       2.000   0.000   0.000  1.00  0.00           C \n"
        "ATOM      4  O   ALA A   1       3.000   0.000   0.000  1.00  0.00           O \n"
        "ATOM      5  N   ALA A   2       4.000   0.000   0.000  1.00  0.00           N \n"
        "ATOM      6  CA  ALA A   2       5.000   0.000   0.000  1.00  0.00           C \n"
        "ATOM      7  N   ALA B   1      10.000   0.000   0.000  1.00  0.00           N \n"
        "ATOM      8  CA  ALA B   1      11.000   0.000   0.000  1.00  0.00           C \n"
        "ATOM      9  C   ALA B   1      12.000   0.000   0.000  1.00  0.00           C \n"
        "ATOM     10  O   ALA B   1      13.000   0.000   0.000  1.00  0.00           O \n"
        "ATOM     11  N   ALA B   2      14.000   0.000   0.000  1.00  0.00           N \n"
        "ATOM     12  CA  ALA B   2      15.000   0.000   0.000  1.00  0.00           C \n";
    test::TempFile input(".pdb", content);
    test::TempFile output(".pdb");

    Molecule molecule(input);
    REQUIRE(molecule.size_body() == 1);
    molecule.save(output);

    auto written = io::detail::pdb::read(output);
    REQUIRE(written.atoms.size() == 12);

    // the two source chains must end up under two different identifiers, in the order they were read
    std::vector<char> chains;
    for (const auto& a : written.atoms) {
        if (chains.empty() || chains.back() != a.chainID) {chains.push_back(a.chainID);}
    }
    REQUIRE(chains.size() == 2);
    CHECK(chains[0] == 'A');
    CHECK(chains[1] == 'B');

    // the residue ids must be preserved, and be unique within each chain
    for (int i = 0; i < 6; ++i) {
        CHECK(written.atoms[i].chainID   == 'A');
        CHECK(written.atoms[i+6].chainID == 'B');
        CHECK(written.atoms[i].resSeq    == written.atoms[i+6].resSeq);
        CHECK(written.atoms[i].resSeq    == (i < 4 ? 1 : 2));
    }

    // the backbone names must survive, since PyMOL and similar tools need them to trace the chains
    std::vector<std::string> expected = {"N", "CA", "C", "O", "N", "CA"};
    for (int i = 0; i < 6; ++i) {
        CHECK(written.atoms[i].name   == expected[i]);
        CHECK(written.atoms[i+6].name == expected[i]);
    }
}

TEST_CASE("PDBWriter: writing multifile pdb") {
    settings::molecule::implicit_hydrogens = false;

    std::vector<AtomFF> atoms(101000);
    for (int i = 0; i < static_cast<int>(atoms.size()); ++i) {
        atoms[i] = AtomFF({1e-4*i, 2e-4*i, 3e-4*i}, form_factor::form_factor_t::C);
    }
    std::vector<Water> waters(100);
    for (int i = 0; i < static_cast<int>(waters.size()); ++i) {
        waters[i] = Water({4e-4*i, 5e-4*i, 6e-4*i});
    }

    Molecule protein({{atoms, waters}});
    io::File base("temp/tests/io/temp_multifile_" + test::detail::unique_tag() + ".pdb");
    protein.save(base);

    io::File part1 = io::File(base).append("_part1");
    io::File part2 = io::File(base).append("_part2");
    REQUIRE(part1.exists());
    REQUIRE(part2.exists());

    // first file
    Molecule protein2(part1);
    REQUIRE(protein2.get_body(0).size_atom() == 100000);

    Molecule protein3(part2);
    REQUIRE(protein3.get_body(0).size_atom() == 1000);
    REQUIRE(protein3.size_water() == 100);

    Molecule protein4(base);
    REQUIRE(protein.size_body() == protein4.size_body());
    for (unsigned int i = 0; i < protein.get_bodies().size(); ++i) {
        if (!protein.get_body(i).equals_content(protein4.get_body(i))) {
            for (unsigned int j = 0; j < protein.get_body(i).size_atom(); ++j) {
                if (!(protein.get_body(i).get_atom(j) == protein4.get_body(i).get_atom(j))) {
                    std::cout << "Difference in body " << i << " atom " << j << std::endl;
                    std::cout << "  Original: " << protein.get_body(i).get_atom(j).coordinates() << " " << protein.get_body(i).get_atom(j).weight() << std::endl;
                    std::cout << "  Loaded:   " << protein4.get_body(i).get_atom(j).coordinates() << " " << protein4.get_body(i).get_atom(j).weight() << std::endl;
                }
            }
        }
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

TEST_CASE("PDBStructure: save") {
    settings::molecule::center = false;
    settings::general::verbose = false;

    auto protein = io::detail::pdb::read("tests/files/2epe.pdb");
    test::TempFile path(".pdb");
    io::detail::pdb::write(protein, path);
    auto protein2 = io::detail::pdb::read(path);
    auto atoms1 = protein.atoms;
    auto atoms2 = protein2.atoms;

    REQUIRE(atoms1.size() == atoms2.size());
    for (unsigned int i = 0; i < atoms1.size(); i++) {
        REQUIRE(atoms1[i].equals_content(atoms2[i]));
    }

    auto waters1 = protein.waters;
    auto waters2 = protein2.waters;
    REQUIRE(waters1.size() == waters2.size());
    for (unsigned int i = 0; i < waters1.size(); i++) {
        REQUIRE(waters1[i].equals_content(waters2[i]));
    }
}