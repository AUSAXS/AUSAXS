#include <catch2/catch_test_macros.hpp>

#include <data/Molecule.h>
#include <data/Body.h>
#include <rigidbody/BodySplitter.h>
#include <io/pdb/PDBStructure.h>
#include <io/pdb/PDBAtom.h>
#include <io/Reader.h>
#include <settings/All.h>

using namespace ausaxs;

TEST_CASE("BodySplitter::split") {
    settings::general::verbose = false;

    auto test_splits = [] (std::string_view file, std::vector<int>&& splits) {
        data::Molecule protein = rigidbody::BodySplitter::split(file, splits);
        io::pdb::PDBStructure data = io::Reader::read(file);
        REQUIRE(protein.size_atom() == data.get_atoms().size());
        
        std::vector<unsigned int> expected_sizes;
        int count = 0;
        for (int split = 0, i = 0; i < static_cast<int>(protein.size_atom()); ++i) {
            auto& ap = data.atoms[i];
            if (ap.resSeq == splits[split]) {
                expected_sizes.push_back(count);
                ++split;
                count = 1;
            } else {
                ++count;
            }
        }
        expected_sizes.push_back(count);

        REQUIRE(protein.size_body() == expected_sizes.size());
        int index = 0;
        for (unsigned int i = 0; i < expected_sizes.size(); ++i) {
            auto& atoms = protein.get_body(i).get_atoms();
            REQUIRE(atoms.size() == expected_sizes[i]);

            if (1 <= i) {
                CHECK(data.atoms[index].resSeq == splits[i-1]);
            }
            if (i < splits.size()) {
                CHECK(data.atoms[index+expected_sizes[i]-1].resSeq == splits[i]-1);
            }
            index += expected_sizes[i];
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