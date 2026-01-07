#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>

#include <io/pdb/PDBAtom.h>
#include <io/pdb/PDBWater.h>
#include <settings/All.h>
#include <constants/Constants.h>

using namespace ausaxs;
using namespace ausaxs::io::pdb;

TEST_CASE("PDBAtom::distance") {
    SECTION("simple") {
        PDBAtom a1({0, 0, 0}, 1, constants::atom_t::N, "GLY", 1);
        PDBAtom a2({1, 0, 0}, 2, constants::atom_t::C, "GLY", 1);
        CHECK(a1.coordinates().distance(a2.coordinates()) == 1);
    }

    SECTION("complex") {
        PDBAtom a1({0, 0, 0}, 1, constants::atom_t::N, "GLY", 1);
        PDBAtom a2({1, 0, 0}, 2, constants::atom_t::C, "GLY", 1);
        PDBAtom a3({1, 1, 0}, 3, constants::atom_t::O, "GLY", 1);
        CHECK(a1.coordinates().distance(a2.coordinates()) == 1);
        CHECK(a1.coordinates().distance(a3.coordinates()) == std::sqrt(2));
    }
}

TEST_CASE("PDBAtom::translate") {
    SECTION("simple") {
        PDBAtom a1({0, 0, 0}, 1, constants::atom_t::N, "GLY", 1);
        a1.coordinates() += Vector3{1, 2, 3};
        CHECK(a1.coordinates() == Vector3{1, 2, 3});
    }

    SECTION("complex") {
        PDBAtom a1({0, 0, 0}, 1, constants::atom_t::N, "GLY", 1);
        a1.coordinates() += Vector3{1, 2, 3};
        CHECK(a1.coordinates() == Vector3<double>{1, 2, 3});
        a1.coordinates() += Vector3{1, 2, 3};
        CHECK(a1.coordinates() == Vector3<double>{2, 4, 6});
    }
}

TEST_CASE("PDBAtom: implicit hydrogens") {
    PDBAtom a(15, "O", "altLoc", "LYS", 'X', 3, "iCode", Vector3<double>{0, 1, 2}, 2.5, 3.5, constants::atom_t::O, "0+");
    CHECK(a.get_mass() == constants::mass::get_mass(constants::atom_t::O) + constants::hydrogen_atoms::residues.get("LYS").get("O", constants::atom_t::O));
    CHECK(a.get_mass() == constants::mass::get_mass(constants::atom_t::O));
}

TEST_CASE("PDBAtom: operators") {
    PDBAtom a1({3, 0, 5}, 2, constants::atom_t::He, "", 3);
    PDBAtom a2 = a1;
    REQUIRE(a1 == a2);
    REQUIRE(a1.equals_content(a2));

    a2 = PDBAtom({0, 4, 1}, 2, constants::atom_t::He, "", 2);
    REQUIRE(a1 != a2);
    REQUIRE(!a1.equals_content(a2));
    REQUIRE(a2 < a1);

    PDBWater w1 = PDBWater({3, 0, 5}, 2, constants::atom_t::He, "", 3);
    PDBWater w2 = w1;
    REQUIRE(w1 == w2);

    w2 = PDBAtom({0, 4, 1}, 2, constants::atom_t::He, "", 2);
    REQUIRE(w1 != w2);
    REQUIRE(w2 < w1);
}

#include <form_factor/FormFactorType.h>
TEST_CASE("PDBAtom: correct_atomic_group_ff") {
    settings::molecule::implicit_hydrogens = true;
    settings::molecule::throw_on_unknown_atom = true;
    PDBAtom atom;

    SECTION("lys") {
        std::string lys1 = "ATOM      1  N   LYS A   1      -3.462  69.119  -8.662  1.00 19.81           N  ";
        std::string lys2 = "ATOM      2  CA  LYS A   1      -2.451  68.681  -9.776  1.00 19.16           C  ";
        std::string lys3 = "ATOM      3  C   LYS A   1      -2.454  67.107  -9.965  1.00 19.10           C  ";
        std::string lys4 = "ATOM      4  O   LYS A   1      -2.418  66.315  -9.018  1.00 16.87           O  ";
        std::string lys5 = "ATOM      5  CB  LYS A   1      -1.010  69.186  -9.464  1.00 21.59           C  ";
        std::string lys6 = "ATOM      6  CG  LYS A   1      -0.034  68.779 -10.377  1.00 25.87           C  ";
        std::string lys7 = "ATOM      7  CD  LYS A   1       1.363  69.238 -10.030  1.00 26.32           C  ";
        std::string lys8 = "ATOM      8  CE  LYS A   1       2.403  68.500 -11.016  1.00 26.04           C  ";
        std::string lys9 = "ATOM      9  NZ  LYS A   1       3.654  69.172 -10.836  1.00 34.18           N  ";
        
        atom.parse_pdb(lys1); atom.add_implicit_hydrogens();
        REQUIRE(atom.atomic_group == constants::atomic_group_t::NH);
        REQUIRE(form_factor::get_type(atom.element, atom.atomic_group) == form_factor::form_factor_t::NH);

        atom.parse_pdb(lys2); atom.add_implicit_hydrogens();
        REQUIRE(atom.atomic_group == constants::atomic_group_t::CH);
        REQUIRE(form_factor::get_type(atom.element, atom.atomic_group) == form_factor::form_factor_t::CH);

        atom.parse_pdb(lys3); atom.add_implicit_hydrogens();
        REQUIRE(atom.atomic_group == constants::atomic_group_t::unknown);
        REQUIRE(form_factor::get_type(atom.element, atom.atomic_group) == form_factor::form_factor_t::C);

        atom.parse_pdb(lys4); atom.add_implicit_hydrogens();
        REQUIRE(atom.atomic_group == constants::atomic_group_t::unknown);
        REQUIRE(form_factor::get_type(atom.element, atom.atomic_group) == form_factor::form_factor_t::O);

        atom.parse_pdb(lys5); atom.add_implicit_hydrogens();
        REQUIRE(atom.atomic_group == constants::atomic_group_t::CH2);
        REQUIRE(form_factor::get_type(atom.element, atom.atomic_group) == form_factor::form_factor_t::CH2);

        atom.parse_pdb(lys6); atom.add_implicit_hydrogens();
        REQUIRE(atom.atomic_group == constants::atomic_group_t::CH2);
        REQUIRE(form_factor::get_type(atom.element, atom.atomic_group) == form_factor::form_factor_t::CH2);

        atom.parse_pdb(lys7); atom.add_implicit_hydrogens();
        REQUIRE(atom.atomic_group == constants::atomic_group_t::CH2);
        REQUIRE(form_factor::get_type(atom.element, atom.atomic_group) == form_factor::form_factor_t::CH2);

        atom.parse_pdb(lys8); atom.add_implicit_hydrogens();
        REQUIRE(atom.atomic_group == constants::atomic_group_t::CH2);
        REQUIRE(form_factor::get_type(atom.element, atom.atomic_group) == form_factor::form_factor_t::CH2);

        atom.parse_pdb(lys9); atom.add_implicit_hydrogens();
        REQUIRE(atom.atomic_group == constants::atomic_group_t::NH3);
        REQUIRE(form_factor::get_type(atom.element, atom.atomic_group) == form_factor::form_factor_t::NH3);
    }

    SECTION("val") {
        std::string val1 = "ATOM     10  N   VAL A   2      -2.619  66.716 -11.199  1.00 19.43           N  ";
        std::string val2 = "ATOM     11  CA  VAL A   2      -2.470  65.345 -11.600  1.00 21.68           C  ";
        std::string val3 = "ATOM     12  C   VAL A   2      -0.988  65.113 -12.076  1.00 21.22           C  ";
        std::string val4 = "ATOM     13  O   VAL A   2      -0.668  65.628 -13.069  1.00 21.74           O  ";
        std::string val5 = "ATOM     14  CB  VAL A   2      -3.483  64.942 -12.686  1.00 19.64           C  ";
        std::string val6 = "ATOM     15  CG1 VAL A   2      -3.247  63.505 -13.005  1.00 17.70           C  ";
        std::string val7 = "ATOM     16  CG2 VAL A   2      -4.940  65.115 -12.243  1.00 19.83           C  ";

        atom.parse_pdb(val1); atom.add_implicit_hydrogens();
        REQUIRE(atom.atomic_group == constants::atomic_group_t::NH);
        REQUIRE(form_factor::get_type(atom.element, atom.atomic_group) == form_factor::form_factor_t::NH);

        atom.parse_pdb(val2); atom.add_implicit_hydrogens();
        REQUIRE(atom.atomic_group == constants::atomic_group_t::CH);
        REQUIRE(form_factor::get_type(atom.element, atom.atomic_group) == form_factor::form_factor_t::CH);

        atom.parse_pdb(val3); atom.add_implicit_hydrogens();
        REQUIRE(atom.atomic_group == constants::atomic_group_t::unknown);
        REQUIRE(form_factor::get_type(atom.element, atom.atomic_group) == form_factor::form_factor_t::C);

        atom.parse_pdb(val4); atom.add_implicit_hydrogens();
        REQUIRE(atom.atomic_group == constants::atomic_group_t::unknown);
        REQUIRE(form_factor::get_type(atom.element, atom.atomic_group) == form_factor::form_factor_t::O);

        atom.parse_pdb(val5); atom.add_implicit_hydrogens();
        REQUIRE(atom.atomic_group == constants::atomic_group_t::CH);
        REQUIRE(form_factor::get_type(atom.element, atom.atomic_group) == form_factor::form_factor_t::CH);

        atom.parse_pdb(val6); atom.add_implicit_hydrogens();
        REQUIRE(atom.atomic_group == constants::atomic_group_t::CH3);
        REQUIRE(form_factor::get_type(atom.element, atom.atomic_group) == form_factor::form_factor_t::CH3);

        atom.parse_pdb(val7); atom.add_implicit_hydrogens();
        REQUIRE(atom.atomic_group == constants::atomic_group_t::CH3);
        REQUIRE(form_factor::get_type(atom.element, atom.atomic_group) == form_factor::form_factor_t::CH3);
    }

    SECTION("phe") {
        std::string phe1 =  "ATOM     17  N   PHE A   3      -0.206  64.328 -11.358  1.00 20.52           N  ";
        std::string phe2 =  "ATOM     18  CA  PHE A   3       1.154  64.049 -11.696  1.00 19.50           C  ";
        std::string phe3 =  "ATOM     19  C   PHE A   3       1.186  63.034 -12.732  1.00 21.23           C  ";
        std::string phe4 =  "ATOM     20  O   PHE A   3       0.286  62.200 -12.856  1.00 22.72           O  ";
        std::string phe5 =  "ATOM     21  CB  PHE A   3       1.929  63.497 -10.445  1.00 19.45           C  ";
        std::string phe6 =  "ATOM     22  CG  PHE A   3       2.500  64.564  -9.596  1.00 19.38           C  ";
        std::string phe7 =  "ATOM     23  CD1 PHE A   3       1.733  65.185  -8.623  1.00 17.20           C  ";
        std::string phe8 =  "ATOM     24  CD2 PHE A   3       3.873  64.910  -9.725  1.00 22.60           C  ";
        std::string phe9 =  "ATOM     25  CE1 PHE A   3       2.290  66.129  -7.768  1.00 21.37           C  ";
        std::string phe10 = "ATOM     26  CE2 PHE A   3       4.425  65.925  -8.883  1.00 26.38           C  ";
        std::string phe11 = "ATOM     27  CZ  PHE A   3       3.575  66.563  -7.911  1.00 24.26           C  ";

        atom.parse_pdb(phe1); atom.add_implicit_hydrogens();
        REQUIRE(atom.atomic_group == constants::atomic_group_t::NH);
        REQUIRE(form_factor::get_type(atom.element, atom.atomic_group) == form_factor::form_factor_t::NH);

        atom.parse_pdb(phe2); atom.add_implicit_hydrogens();
        REQUIRE(atom.atomic_group == constants::atomic_group_t::CH);
        REQUIRE(form_factor::get_type(atom.element, atom.atomic_group) == form_factor::form_factor_t::CH);

        atom.parse_pdb(phe3); atom.add_implicit_hydrogens();
        REQUIRE(atom.atomic_group == constants::atomic_group_t::unknown);
        REQUIRE(form_factor::get_type(atom.element, atom.atomic_group) == form_factor::form_factor_t::C);

        atom.parse_pdb(phe4); atom.add_implicit_hydrogens();
        REQUIRE(atom.atomic_group == constants::atomic_group_t::unknown);
        REQUIRE(form_factor::get_type(atom.element, atom.atomic_group) == form_factor::form_factor_t::O);

        atom.parse_pdb(phe5); atom.add_implicit_hydrogens();
        REQUIRE(atom.atomic_group == constants::atomic_group_t::CH2);
        REQUIRE(form_factor::get_type(atom.element, atom.atomic_group) == form_factor::form_factor_t::CH2);

        atom.parse_pdb(phe6); atom.add_implicit_hydrogens();
        REQUIRE(atom.atomic_group == constants::atomic_group_t::unknown);
        REQUIRE(form_factor::get_type(atom.element, atom.atomic_group) == form_factor::form_factor_t::C);

        atom.parse_pdb(phe7); atom.add_implicit_hydrogens();
        REQUIRE(atom.atomic_group == constants::atomic_group_t::CH);
        REQUIRE(form_factor::get_type(atom.element, atom.atomic_group) == form_factor::form_factor_t::CH);

        atom.parse_pdb(phe8); atom.add_implicit_hydrogens();
        REQUIRE(atom.atomic_group == constants::atomic_group_t::CH);
        REQUIRE(form_factor::get_type(atom.element, atom.atomic_group) == form_factor::form_factor_t::CH);

        atom.parse_pdb(phe9); atom.add_implicit_hydrogens();
        REQUIRE(atom.atomic_group == constants::atomic_group_t::CH);
        REQUIRE(form_factor::get_type(atom.element, atom.atomic_group) == form_factor::form_factor_t::CH);

        atom.parse_pdb(phe10); atom.add_implicit_hydrogens();
        REQUIRE(atom.atomic_group == constants::atomic_group_t::CH);
        REQUIRE(form_factor::get_type(atom.element, atom.atomic_group) == form_factor::form_factor_t::CH);

        atom.parse_pdb(phe11); atom.add_implicit_hydrogens();
        REQUIRE(atom.atomic_group == constants::atomic_group_t::CH);
        REQUIRE(form_factor::get_type(atom.element, atom.atomic_group) == form_factor::form_factor_t::CH);
    }

    SECTION("gly") {
        std::string gly1 = "ATOM     28  N   GLY A   4       2.287  63.055 -13.488  1.00 21.95           N  ";
        std::string gly2 = "ATOM     29  CA  GLY A   4       2.605  61.971 -14.393  1.00 19.79           C  ";
        std::string gly3 = "ATOM     30  C   GLY A   4       3.475  60.975 -13.566  1.00 19.47           C  ";
        std::string gly4 = "ATOM     31  O   GLY A   4       3.990  61.318 -12.551  1.00 16.69           O  ";

        atom.parse_pdb(gly1); atom.add_implicit_hydrogens();
        REQUIRE(atom.atomic_group == constants::atomic_group_t::NH);
        REQUIRE(form_factor::get_type(atom.element, atom.atomic_group) == form_factor::form_factor_t::NH);

        atom.parse_pdb(gly2); atom.add_implicit_hydrogens();
        REQUIRE(atom.atomic_group == constants::atomic_group_t::CH2);
        REQUIRE(form_factor::get_type(atom.element, atom.atomic_group) == form_factor::form_factor_t::CH2);

        atom.parse_pdb(gly3); atom.add_implicit_hydrogens();
        REQUIRE(atom.atomic_group == constants::atomic_group_t::unknown);
        REQUIRE(form_factor::get_type(atom.element, atom.atomic_group) == form_factor::form_factor_t::C);

        atom.parse_pdb(gly4); atom.add_implicit_hydrogens();
        REQUIRE(atom.atomic_group == constants::atomic_group_t::unknown);
        REQUIRE(form_factor::get_type(atom.element, atom.atomic_group) == form_factor::form_factor_t::O);
    }

    SECTION("met") {
        std::string met1 = "ATOM     42  N   MET A   6       2.683  -9.695  -4.055  1.00 35.86           N  ";
        std::string met2 = "ATOM     43  CA  MET A   6       2.271 -11.076  -4.245  1.00 38.24           C  ";
        std::string met3 = "ATOM     44  C   MET A   6       3.262 -12.007  -3.567  1.00 35.64           C  ";
        std::string met4 = "ATOM     45  O   MET A   6       4.477 -11.842  -3.708  1.00 35.28           O  ";
        std::string met5 = "ATOM     46  CB  MET A   6       2.177 -11.397  -5.740  1.00 47.86           C  ";
        std::string met6 = "ATOM     47  CG  MET A   6       1.540 -12.723  -6.078  1.00 54.30           C  ";
        std::string met7 = "ATOM     48  SD  MET A   6       1.467 -12.932  -7.867  1.00 55.60           S  ";
        std::string met8 = "ATOM     49  CE  MET A   6       0.762 -11.361  -8.347  1.00 49.93           C  ";

        atom.parse_pdb(met1); atom.add_implicit_hydrogens();
        REQUIRE(atom.atomic_group == constants::atomic_group_t::NH);
        REQUIRE(form_factor::get_type(atom.element, atom.atomic_group) == form_factor::form_factor_t::NH);

        atom.parse_pdb(met2); atom.add_implicit_hydrogens();
        REQUIRE(atom.atomic_group == constants::atomic_group_t::CH);
        REQUIRE(form_factor::get_type(atom.element, atom.atomic_group) == form_factor::form_factor_t::CH);

        atom.parse_pdb(met3); atom.add_implicit_hydrogens();
        REQUIRE(atom.atomic_group == constants::atomic_group_t::unknown);
        REQUIRE(form_factor::get_type(atom.element, atom.atomic_group) == form_factor::form_factor_t::C);

        atom.parse_pdb(met4); atom.add_implicit_hydrogens();
        REQUIRE(atom.atomic_group == constants::atomic_group_t::unknown);
        REQUIRE(form_factor::get_type(atom.element, atom.atomic_group) == form_factor::form_factor_t::O);

        atom.parse_pdb(met5); atom.add_implicit_hydrogens();
        REQUIRE(atom.atomic_group == constants::atomic_group_t::CH2);
        REQUIRE(form_factor::get_type(atom.element, atom.atomic_group) == form_factor::form_factor_t::CH2);

        atom.parse_pdb(met6); atom.add_implicit_hydrogens();
        REQUIRE(atom.atomic_group == constants::atomic_group_t::CH2);
        REQUIRE(form_factor::get_type(atom.element, atom.atomic_group) == form_factor::form_factor_t::CH2);

        atom.parse_pdb(met7); atom.add_implicit_hydrogens();
        REQUIRE(atom.atomic_group == constants::atomic_group_t::unknown);
        REQUIRE(form_factor::get_type(atom.element, atom.atomic_group) == form_factor::form_factor_t::S);

        atom.parse_pdb(met8); atom.add_implicit_hydrogens();
        REQUIRE(atom.atomic_group == constants::atomic_group_t::CH3);
        REQUIRE(form_factor::get_type(atom.element, atom.atomic_group) == form_factor::form_factor_t::CH3);
    }
}