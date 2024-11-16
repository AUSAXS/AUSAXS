#include <catch2/catch_test_macros.hpp>

#include <data/Molecule.h>
#include <data/Body.h>
#include <data/record/Atom.h>
#include <rigidbody/BodySplitter.h>
#include <settings/All.h>

using namespace ausaxs;

TEST_CASE("BodySplitter::split") {
    settings::general::verbose = false;

    auto test_splits = [] (std::string_view file, std::vector<int>&& splits) {
        data::Molecule protein = rigidbody::BodySplitter::split(file, splits);
        
        std::vector<unsigned int> expected_sizes;
        int count = 0;
        for (int split = 0; auto& a : protein.get_atoms()) {
            if (a.get_residue_sequence_number() == splits[split]) {
                expected_sizes.push_back(count);
                ++split;
                count = 1;
            } else {
                ++count;
            }
        }
        expected_sizes.push_back(count);

        REQUIRE(protein.size_body() == expected_sizes.size());
        for (unsigned int i = 0; i < expected_sizes.size(); ++i) {
            auto& atoms = protein.get_body(i).get_atoms();
            REQUIRE(atoms.size() == expected_sizes[i]);

            if (1 <= i) {
                CHECK(atoms.front().get_residue_sequence_number() == splits[i-1]);
            }
            if (i < splits.size()) {
                CHECK(atoms.back().get_residue_sequence_number() == splits[i]-1);
            }
        }
    };

    SECTION("LAR1-2") {
        SECTION("1 split") {
            test_splits("tests/files/LAR1-2.pdb", {9});
        }

        SECTION("2 splits") {
            test_splits("tests/files/LAR1-2.pdb", {9, 99});
        }

        SECTION("many splits") {
            test_splits("tests/files/LAR1-2.pdb", {5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100, 200});
        }
    }

    SECTION("2epe") {
        SECTION("1 split") {
            test_splits("tests/files/2epe.pdb", {10});
        }

        SECTION("2 splits") {
            test_splits("tests/files/2epe.pdb", {10, 50});
        }

        SECTION("many splits") {
            test_splits("tests/files/2epe.pdb", {5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100});
        }
    }

    SECTION("SASDJG5") {
        SECTION("1 split") {
            test_splits("tests/files/SASDJG5.pdb", {10});
        }

        SECTION("2 splits") {
            test_splits("tests/files/SASDJG5.pdb", {10, 50});
        }

        SECTION("many splits") {
            test_splits("tests/files/SASDJG5.pdb", {5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100});
        }

        SECTION("by chain") {
            data::Molecule protein = rigidbody::BodySplitter::split("tests/files/SASDJG5.pdb");
            REQUIRE(protein.size_body() == 2);
            REQUIRE(protein.get_body(0).size_atom() == 2367);
            REQUIRE(protein.get_body(1).size_atom() == 2367);
        }
    }
}