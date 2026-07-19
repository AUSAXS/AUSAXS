// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <catch2/catch_test_macros.hpp>

#include <rigidbody/sequencer/detail/SequenceParser.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/Rigidbody.h>
#include <rigidbody/detail/SystemSpecification.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <settings/All.h>
#include <io/ExistingFile.h>

#include <fstream>

using namespace ausaxs;
using namespace ausaxs::rigidbody;
using namespace ausaxs::rigidbody::sequencer;

struct SequenceParserMergeFixture {
    SequenceParserMergeFixture() {
        settings::general::verbose = false;
        settings::molecule::implicit_hydrogens = false;
        settings::grid::min_bins = 250;
        settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::None;
    }

    std::unique_ptr<Sequencer> parse(const std::string& content) {
        static int counter = 0;
        std::string path = "/tmp/ausaxs_seq_merge_test_" + std::to_string(counter++) + ".conf";
        std::ofstream f(path);
        f << content;
        f.close();
        SequenceParser parser;
        return parser.parse_file(path);
    }
};

TEST_CASE_METHOD(SequenceParserMergeFixture, "SequenceParser::MergeElement") {
    SECTION("merging two bodies combines their atoms into one body") {
        auto baseline = parse(
            "load {\n"
            "    pdb tests/files/SASDJG5_single.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
        );
        REQUIRE(baseline != nullptr);
        auto single_body_atoms = baseline->_get_rigidbody()->molecule.get_body(0).size_atom();

        auto seq = parse(
            "load {\n"
            "    pdb tests/files/SASDJG5_single.pdb tests/files/SASDJG5_single.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
            "merge b1 b2\n"
        );
        REQUIRE(seq != nullptr);
        auto rb = seq->_get_rigidbody();
        REQUIRE(rb != nullptr);

        REQUIRE(rb->molecule.size_body() == 1);
        CHECK(rb->molecule.get_body(0).size_atom() == 2*single_body_atoms);
    }

    SECTION("merging more than one body into the first works in a single statement") {
        auto seq = parse(
            "load {\n"
            "    pdb tests/files/SASDJG5_single.pdb tests/files/SASDJG5_single.pdb tests/files/SASDJG5_single.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
            "merge b1 b2 b3\n"
        );
        REQUIRE(seq != nullptr);
        auto rb = seq->_get_rigidbody();
        REQUIRE(rb != nullptr);
        CHECK(rb->molecule.size_body() == 1);
    }

    SECTION("names of merged-away bodies are removed from the setup") {
        auto seq = parse(
            "load {\n"
            "    pdb tests/files/SASDJG5_single.pdb tests/files/SASDJG5_single.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
            "merge b1 b2\n"
        );
        REQUIRE(seq != nullptr);
        const auto& names = seq->setup()._get_body_names();
        CHECK(names.contains("b1"));
        CHECK_FALSE(names.contains("b2"));
    }

    SECTION("surviving body names are reindexed after bodies before them are merged away") {
        auto seq = parse(
            "load {\n"
            "    pdb tests/files/SASDJG5_single.pdb tests/files/SASDJG5_single.pdb tests/files/SASDJG5_single.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
            "merge b1 b2\n"
            "symmetry {\n"
            "    b3 c2\n"
            "}\n"
        );
        REQUIRE(seq != nullptr);
        auto rb = seq->_get_rigidbody();
        REQUIRE(rb != nullptr);
        REQUIRE(rb->molecule.size_body() == 2);

        // b3 was originally at index 2; after b2 (index 1) is merged away it should have shifted down to index 1
        CHECK(rb->molecule.get_body(0).size_symmetry() == 0);
        CHECK(rb->molecule.get_body(1).size_symmetry() == 1);
    }

    SECTION("merging a body into itself is rejected") {
        CHECK_THROWS(parse(
            "load {\n"
            "    pdb tests/files/SASDJG5_single.pdb tests/files/SASDJG5_single.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
            "merge b1 b1\n"
        ));
    }

    SECTION("merging an unknown body name is rejected") {
        CHECK_THROWS(parse(
            "load {\n"
            "    pdb tests/files/SASDJG5_single.pdb tests/files/SASDJG5_single.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
            "merge b1 doesnotexist\n"
        ));
    }

    SECTION("a bare merge with only one argument is rejected") {
        CHECK_THROWS(parse(
            "load {\n"
            "    pdb tests/files/SASDJG5_single.pdb tests/files/SASDJG5_single.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
            "merge b1\n"
        ));
    }
}
