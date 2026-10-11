#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <data/Body.h>
#include <data/Molecule.h>
#include <io/detail/structure/PDBReader.h>
#include <io/detail/structure/PDBWriter.h>
#include <io/pdb/PDBStructure.h>
#include <settings/All.h>

#include <support/temp_file.h>

#include <cmath>
#include <string>
#include <vector>

using namespace ausaxs;
using namespace ausaxs::data;
using namespace ausaxs::io::pdb;

TEST_CASE("PDBStructure: save") {
    settings::general::verbose = false;

    auto protein = io::detail::pdb::read("tests/files/2epe.pdb");
    test::TempFile path(".pdb");
    io::detail::pdb::write(protein, path);
    auto protein2 = io::detail::pdb::read(path);
    auto atoms1 = protein.atoms;
    auto atoms2 = protein2.atoms;

    REQUIRE(atoms1.size() == atoms2.size());
    for (int i = 0; i < static_cast<int>(atoms1.size()); i++) {
        REQUIRE(atoms1[i].equals_content(atoms2[i]));
    }

    auto waters1 = protein.waters;
    auto waters2 = protein2.waters;
    REQUIRE(waters1.size() == waters2.size());
    for (int i = 0; i < static_cast<int>(waters1.size()); i++) {
        REQUIRE(waters1[i].equals_content(waters2[i]));
    }
}

TEST_CASE("PDBStructure: a written structure can be read back unchanged") {
    settings::general::verbose = false;
    settings::molecule::use_occupancy = true;
    settings::molecule::implicit_hydrogens = true;
    settings::molecule::allow_unknown_residues = false;

    // The reduced representation keeps only coordinates, form factors and charges, so everything the residue table is keyed on - the residue name, the atom
    // name and the occupancy - has to come back out of the metadata when the structure is written. Without it every residue is written as "UNK", which makes
    // the implicit hydrogen assignment throw on the next read, and forcing it through instead drops the implicit hydrogens of every atom.
    SECTION("hand-written structure") {
        std::string content =
            "ATOM      1  N   LYS A   1       0.000   0.000   0.000  1.00  0.00           N \n"
            "ATOM      2  CA  LYS A   1       1.000   0.000   0.000  1.00  0.00           C \n"
            "ATOM      3  CB  LYS A   1       2.000   0.000   0.000  1.00  0.00           C \n"
            "ATOM      4  CG  LYS A   1       3.000   0.000   0.000  0.50  0.00           C \n"
            "ATOM      5  C   LYS A   1       4.000   0.000   0.000  1.00  0.00           C \n"
            "ATOM      6  O   LYS A   1       5.000   0.000   0.000  1.00  0.00           O \n";
        test::TempFile input(".pdb", content);
        test::TempFile output(".pdb");

        Molecule molecule(input);
        molecule.save(output);

        auto written = io::detail::pdb::read(output);
        REQUIRE(written.atoms.size() == 6);

        std::vector<std::string> names = {"N", "CA", "CB", "CG", "C", "O"};
        for (int i = 0; i < 6; ++i) {
            CHECK(written.atoms[i].resName == "LYS");
            CHECK(written.atoms[i].name    == names[i]);
        }

        // the partially occupied atom contributes half its charge; writing 1.00 for it would inflate that on the next read
        CHECK_THAT(written.atoms[3].occupancy, Catch::Matchers::WithinAbs(0.5, 1e-6));

        Molecule reloaded(output);
        REQUIRE(reloaded.size_atom() == molecule.size_atom());
        for (int i = 0; i < reloaded.size_atom(); ++i) {
            CHECK(reloaded.get_body(0).get_atom(i).form_factor_type() == molecule.get_body(0).get_atom(i).form_factor_type());
        }
        CHECK_THAT(reloaded.get_total_atomic_charge(), Catch::Matchers::WithinRel(molecule.get_total_atomic_charge(), 1e-12));
    }

    SECTION("real structure") {
        test::TempFile output(".pdb");
        Molecule molecule("tests/files/2epe.pdb");
        molecule.save(output);

        Molecule reloaded(output);
        REQUIRE(reloaded.size_atom() == molecule.size_atom());
        for (int i = 0; i < reloaded.size_atom(); ++i) {
            CHECK(reloaded.get_body(0).get_atom(i).form_factor_type() == molecule.get_body(0).get_atom(i).form_factor_type());
        }
        CHECK_THAT(reloaded.get_total_atomic_charge(), Catch::Matchers::WithinRel(molecule.get_total_atomic_charge(), 1e-12));
        CHECK_THAT(reloaded.get_absolute_mass(),       Catch::Matchers::WithinRel(molecule.get_absolute_mass(),       1e-12));
    }

    SECTION("elements outside the tabulated form factors") {
        settings::molecule::implicit_hydrogens = false;
        std::string content =
            "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N \n"
            "ATOM      2  CA  ALA A   1       1.000   0.000   0.000  1.00  0.00           C \n"
            "HETATM    3  P   PO4 A 101       4.000   0.000   0.000  1.00  0.00           P \n"
            "HETATM    4 ZN    ZN A 102       7.000   0.000   0.000  1.00  0.00          ZN \n"
            "HETATM    5 CL    CL A 103      10.000   0.000   0.000  0.50  0.00          CL \n";
        test::TempFile input(".pdb", content);
        test::TempFile output(".pdb");

        Molecule molecule(input);
        molecule.save(output);

        auto written = io::detail::pdb::read(output);
        REQUIRE(written.atoms.size() == 5);
        CHECK(written.atoms[2].element == constants::atom_t::P);
        CHECK(written.atoms[3].element == constants::atom_t::Zn);
        CHECK(written.atoms[4].element == constants::atom_t::Cl);

        Molecule reloaded(output);
        REQUIRE(reloaded.size_atom() == molecule.size_atom());
        for (int i = 0; i < reloaded.size_atom(); ++i) {
            CHECK(reloaded.get_body(0).get_atom(i).form_factor_type() == molecule.get_body(0).get_atom(i).form_factor_type());
            CHECK_THAT(reloaded.get_body(0).get_atom(i).weight(), Catch::Matchers::WithinRel(molecule.get_body(0).get_atom(i).weight(), 1e-6));
        }
        settings::molecule::implicit_hydrogens = true;
    }
}

namespace {
    // reduce a structure consisting of a single atom parsed from a PDB line
    AtomFF reduce_single(const std::string& line) {
        PDBAtom atom;
        atom.parse_pdb(line);
        return PDBStructure({atom}, {}).reduced_representation().atoms[0];
    }
}

TEST_CASE("PDBStructure: correct_atomic_group_ff") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = true;
    settings::molecule::allow_unknown_atoms = false;
    auto ff_type = [] (const std::string& line) {return reduce_single(line).form_factor_type();};

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
        
        REQUIRE(ff_type(lys1) == form_factor::form_factor_t::NH);
        REQUIRE(ff_type(lys2) == form_factor::form_factor_t::CH);
        REQUIRE(ff_type(lys3) == form_factor::form_factor_t::C);
        REQUIRE(ff_type(lys4) == form_factor::form_factor_t::O);
        REQUIRE(ff_type(lys5) == form_factor::form_factor_t::CH2);
        REQUIRE(ff_type(lys6) == form_factor::form_factor_t::CH2);
        REQUIRE(ff_type(lys7) == form_factor::form_factor_t::CH2);
        REQUIRE(ff_type(lys8) == form_factor::form_factor_t::CH2);
        REQUIRE(ff_type(lys9) == form_factor::form_factor_t::NH3);
    }

    SECTION("val") {
        std::string val1 = "ATOM     10  N   VAL A   2      -2.619  66.716 -11.199  1.00 19.43           N  ";
        std::string val2 = "ATOM     11  CA  VAL A   2      -2.470  65.345 -11.600  1.00 21.68           C  ";
        std::string val3 = "ATOM     12  C   VAL A   2      -0.988  65.113 -12.076  1.00 21.22           C  ";
        std::string val4 = "ATOM     13  O   VAL A   2      -0.668  65.628 -13.069  1.00 21.74           O  ";
        std::string val5 = "ATOM     14  CB  VAL A   2      -3.483  64.942 -12.686  1.00 19.64           C  ";
        std::string val6 = "ATOM     15  CG1 VAL A   2      -3.247  63.505 -13.005  1.00 17.70           C  ";
        std::string val7 = "ATOM     16  CG2 VAL A   2      -4.940  65.115 -12.243  1.00 19.83           C  ";

        REQUIRE(ff_type(val1) == form_factor::form_factor_t::NH);
        REQUIRE(ff_type(val2) == form_factor::form_factor_t::CH);
        REQUIRE(ff_type(val3) == form_factor::form_factor_t::C);
        REQUIRE(ff_type(val4) == form_factor::form_factor_t::O);
        REQUIRE(ff_type(val5) == form_factor::form_factor_t::CH);
        REQUIRE(ff_type(val6) == form_factor::form_factor_t::CH3);
        REQUIRE(ff_type(val7) == form_factor::form_factor_t::CH3);
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

        REQUIRE(ff_type(phe1) == form_factor::form_factor_t::NH);
        REQUIRE(ff_type(phe2) == form_factor::form_factor_t::CH);
        REQUIRE(ff_type(phe3) == form_factor::form_factor_t::C);
        REQUIRE(ff_type(phe4) == form_factor::form_factor_t::O);
        REQUIRE(ff_type(phe5) == form_factor::form_factor_t::CH2);
        REQUIRE(ff_type(phe6) == form_factor::form_factor_t::C);
        REQUIRE(ff_type(phe7) == form_factor::form_factor_t::CH);
        REQUIRE(ff_type(phe8) == form_factor::form_factor_t::CH);
        REQUIRE(ff_type(phe9) == form_factor::form_factor_t::CH);
        REQUIRE(ff_type(phe10) == form_factor::form_factor_t::CH);
        REQUIRE(ff_type(phe11) == form_factor::form_factor_t::CH);
    }

    SECTION("gly") {
        std::string gly1 = "ATOM     28  N   GLY A   4       2.287  63.055 -13.488  1.00 21.95           N  ";
        std::string gly2 = "ATOM     29  CA  GLY A   4       2.605  61.971 -14.393  1.00 19.79           C  ";
        std::string gly3 = "ATOM     30  C   GLY A   4       3.475  60.975 -13.566  1.00 19.47           C  ";
        std::string gly4 = "ATOM     31  O   GLY A   4       3.990  61.318 -12.551  1.00 16.69           O  ";

        REQUIRE(ff_type(gly1) == form_factor::form_factor_t::NH);
        REQUIRE(ff_type(gly2) == form_factor::form_factor_t::CH2);
        REQUIRE(ff_type(gly3) == form_factor::form_factor_t::C);
        REQUIRE(ff_type(gly4) == form_factor::form_factor_t::O);
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

        REQUIRE(ff_type(met1) == form_factor::form_factor_t::NH);
        REQUIRE(ff_type(met2) == form_factor::form_factor_t::CH);
        REQUIRE(ff_type(met3) == form_factor::form_factor_t::C);
        REQUIRE(ff_type(met4) == form_factor::form_factor_t::O);
        REQUIRE(ff_type(met5) == form_factor::form_factor_t::CH2);
        REQUIRE(ff_type(met6) == form_factor::form_factor_t::CH2);
        REQUIRE(ff_type(met7) == form_factor::form_factor_t::S);
        REQUIRE(ff_type(met8) == form_factor::form_factor_t::CH3);
    }
}

TEST_CASE("PDBStructure: implicit hydrogens are counted once") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = true;
    settings::molecule::allow_unknown_atoms = false;

    // a group form factor already carries its hydrogens, so its charge is the whole group: CA of LYS is CH, 7 electrons
    auto atom = reduce_single("ATOM      2  CA  LYS A   1      -2.451  68.681  -9.776  1.00 19.16           C  ");
    REQUIRE(atom.form_factor_type() == form_factor::form_factor_t::CH);
    CHECK(atom.weight() == constants::charge::get_ff_charge(form_factor::form_factor_t::CH));
    CHECK(std::round(atom.weight()) == 7);

    // likewise SG of CYS is SH, 17 electrons
    atom = reduce_single("ATOM      6  SG  CYS A   1       0.000   0.000   0.000  1.00 19.16           S  ");
    REQUIRE(atom.form_factor_type() == form_factor::form_factor_t::SH);
    CHECK(std::round(atom.weight()) == 17);
}
