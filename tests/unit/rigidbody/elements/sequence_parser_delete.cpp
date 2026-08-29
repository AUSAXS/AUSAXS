// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <catch2/catch_test_macros.hpp>

#include <rigidbody/sequencer/detail/SequenceParser.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/Rigidbody.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <settings/All.h>
#include <io/ExistingFile.h>

#include <support/temp_file.h>


using namespace ausaxs;
using namespace ausaxs::rigidbody;
using namespace ausaxs::rigidbody::sequencer;

struct SequenceParserDeleteFixture {
    SequenceParserDeleteFixture() {
        settings::general::verbose = false;
        settings::molecule::implicit_hydrogens = false;
        settings::grid::min_bins = 250;
        settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::None;
    }

    std::unique_ptr<Sequencer> parse(const std::string& content) {
        test::TempFile config(".conf", content);
        SequenceParser parser;
        return parser.parse_file(config);
    }
};

TEST_CASE_METHOD(SequenceParserDeleteFixture, "SequenceParser::DeleteElement") {
    SECTION("deleting a body removes it from the molecule") {
        auto seq = parse(
            "load {\n"
            "    pdb tests/files/SASDJG5_single.pdb tests/files/SASDJG5_single.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
            "delete b2\n"
        );
        REQUIRE(seq != nullptr);
        auto rb = seq->_get_rigidbody();
        REQUIRE(rb != nullptr);
        CHECK(rb->molecule.size_body() == 1);
    }

    SECTION("multiple bodies can be deleted in a single statement") {
        auto seq = parse(
            "load {\n"
            "    pdb tests/files/SASDJG5_single.pdb tests/files/SASDJG5_single.pdb tests/files/SASDJG5_single.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
            "delete b1 b2\n"
        );
        REQUIRE(seq != nullptr);
        auto rb = seq->_get_rigidbody();
        REQUIRE(rb != nullptr);
        CHECK(rb->molecule.size_body() == 1);
    }

    SECTION("deleted body names are removed from the setup") {
        auto seq = parse(
            "load {\n"
            "    pdb tests/files/SASDJG5_single.pdb tests/files/SASDJG5_single.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
            "delete b2\n"
        );
        REQUIRE(seq != nullptr);
        const auto& names = seq->setup()._body_name_registry();
        CHECK(names.contains("b1"));
        CHECK_FALSE(names.contains("b2"));
    }

    SECTION("surviving body names are reindexed after an earlier body is deleted") {
        auto seq = parse(
            "load {\n"
            "    pdb tests/files/SASDJG5_single.pdb tests/files/SASDJG5_single.pdb tests/files/SASDJG5_single.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
            "delete b1\n"
            "symmetry b3 c2\n"
        );
        REQUIRE(seq != nullptr);
        auto rb = seq->_get_rigidbody();
        REQUIRE(rb != nullptr);
        REQUIRE(rb->molecule.size_body() == 2);

        // b3 was originally at index 2; after b1 (index 0) is deleted it should have shifted down to index 1
        CHECK(rb->molecule.get_body(0).size_symmetry() == 0);
        CHECK(rb->molecule.get_body(1).size_symmetry() == 1);
    }

    SECTION("list separators like commas are ignored in list arguments") {
        auto seq = parse(
            "load {\n"
            "    pdb tests/files/SASDJG5_single.pdb tests/files/SASDJG5_single.pdb tests/files/SASDJG5_single.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
            "delete b1, b2\n"
        );
        REQUIRE(seq != nullptr);
        auto rb = seq->_get_rigidbody();
        REQUIRE(rb != nullptr);
        CHECK(rb->molecule.size_body() == 1);
    }

    SECTION("deleting an unknown body name is rejected") {
        CHECK_THROWS(parse(
            "load {\n"
            "    pdb tests/files/SASDJG5_single.pdb tests/files/SASDJG5_single.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
            "delete doesnotexist\n"
        ));
    }

    SECTION("deleting the same body twice in one statement is rejected") {
        CHECK_THROWS(parse(
            "load {\n"
            "    pdb tests/files/SASDJG5_single.pdb tests/files/SASDJG5_single.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
            "delete b1 b1\n"
        ));
    }

    SECTION("deleting every body is rejected") {
        CHECK_THROWS(parse(
            "load {\n"
            "    pdb tests/files/SASDJG5_single.pdb tests/files/SASDJG5_single.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
            "delete b1 b2\n"
        ));
    }

    SECTION("a bare delete with no arguments is rejected") {
        CHECK_THROWS(parse(
            "load {\n"
            "    pdb tests/files/SASDJG5_single.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
            "delete\n"
        ));
    }
}
