#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <io/detail/structure/PDBReader.h>
#include <io/detail/structure/PDBWriter.h>
#include <io/pdb/PDBStructure.h>
#include <settings/All.h>

#include <support/temp_file.h>

#include <string>
#include <vector>

using namespace ausaxs;
using namespace ausaxs::data;
using namespace ausaxs::io::pdb;

TEST_CASE("PDBReader::read") {
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

        CHECK_THAT(atoms[0].effective_charge, Catch::Matchers::WithinRel(constants::charge::get_ff_charge(atoms[0].get_form_factor_type()), 1e-12));
        CHECK(atoms[0].atomic_group == constants::atomic_group_t::NH);

        CHECK_THAT(atoms[1].effective_charge, Catch::Matchers::WithinRel(constants::charge::get_ff_charge(atoms[1].get_form_factor_type()), 1e-12));
        CHECK(atoms[1].atomic_group == constants::atomic_group_t::CH);

        CHECK_THAT(atoms[2].effective_charge, Catch::Matchers::WithinRel(constants::charge::get_ff_charge(atoms[2].get_form_factor_type()) + 0, 1e-12));
        CHECK(atoms[2].atomic_group == constants::atomic_group_t::unknown);

        CHECK_THAT(atoms[3].effective_charge, Catch::Matchers::WithinRel(constants::charge::get_ff_charge(atoms[3].get_form_factor_type()) + 0, 1e-12));
        CHECK(atoms[3].atomic_group == constants::atomic_group_t::unknown);

        CHECK_THAT(atoms[4].effective_charge, Catch::Matchers::WithinRel(constants::charge::get_ff_charge(atoms[4].get_form_factor_type()), 1e-12));
        CHECK(atoms[4].atomic_group == constants::atomic_group_t::CH2);

        CHECK_THAT(atoms[5].effective_charge, Catch::Matchers::WithinRel(constants::charge::get_ff_charge(atoms[5].get_form_factor_type()), 1e-12));
        CHECK(atoms[5].atomic_group == constants::atomic_group_t::CH2);

        CHECK_THAT(atoms[6].effective_charge, Catch::Matchers::WithinRel(constants::charge::get_ff_charge(atoms[6].get_form_factor_type()), 1e-12));
        CHECK(atoms[6].atomic_group == constants::atomic_group_t::CH2);

        CHECK_THAT(atoms[7].effective_charge, Catch::Matchers::WithinRel(constants::charge::get_ff_charge(atoms[7].get_form_factor_type()), 1e-12));
        CHECK(atoms[7].atomic_group == constants::atomic_group_t::CH2);

        CHECK_THAT(atoms[8].effective_charge, Catch::Matchers::WithinRel(constants::charge::get_ff_charge(atoms[8].get_form_factor_type()), 1e-12));
        CHECK(atoms[8].atomic_group == constants::atomic_group_t::NH3);
    }

    SECTION("disabled") {
        auto protein = generate_molecule();

        for (auto a : protein.atoms) {
            CHECK_THAT(a.effective_charge, Catch::Matchers::WithinRel(constants::charge::get_ff_charge(a.get_form_factor_type()), 1e-12));
            CHECK(a.atomic_group == constants::atomic_group_t::unknown);
        }
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
