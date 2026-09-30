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

TEST_CASE("PDBWriter: multi-chain structures keep their chains apart") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;

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
    for (int i = 0; i < static_cast<int>(protein.get_bodies().size()); ++i) {
        if (!protein.get_body(i).equals_content(protein4.get_body(i))) {
            for (int j = 0; j < protein.get_body(i).size_atom(); ++j) {
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
