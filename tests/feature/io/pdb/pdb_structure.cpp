#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <data/Body.h>
#include <data/Molecule.h>
#include <io/detail/structure/PDBReader.h>
#include <io/detail/structure/PDBWriter.h>
#include <io/pdb/PDBStructure.h>
#include <settings/All.h>

#include <support/temp_file.h>

#include <iostream>
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
    // add_implicit_hydrogens throw on the next read, and forcing it through instead drops the implicit hydrogens of every atom.
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
}
