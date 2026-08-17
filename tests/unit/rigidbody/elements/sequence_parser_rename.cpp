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

struct SequenceParserRenameFixture {
    SequenceParserRenameFixture() {
        settings::general::verbose = false;
        settings::molecule::implicit_hydrogens = false;
        settings::grid::min_bins = 250;
        settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::None;
    }

    std::unique_ptr<Sequencer> parse(const std::string& content) {
        test::TempFile config("ausaxs_seq_rename_test", ".conf", content);
        SequenceParser parser;
        return parser.parse_file(config);
    }
};

TEST_CASE_METHOD(SequenceParserRenameFixture, "SequenceParser::RenameElement") {
    SECTION("renaming a body updates the setup's name map") {
        auto seq = parse(
            "load {\n"
            "    pdb tests/files/SASDJG5_single.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
            "rename b1 core\n"
        );
        REQUIRE(seq != nullptr);
        const auto& names = seq->setup()._body_name_registry();
        CHECK(names.contains("b1")); // the default name is permanent and is never removed by a rename
        CHECK(names.contains("core"));
    }

    SECTION("renaming again replaces the previous alias instead of accumulating it") {
        auto seq = parse(
            "load {\n"
            "    pdb tests/files/SASDJG5_single.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
            "rename b1 core\n"
            "rename core shell\n"
        );
        REQUIRE(seq != nullptr);
        const auto& names = seq->setup()._body_name_registry();
        CHECK(names.contains("b1"));
        CHECK_FALSE(names.contains("core"));
        CHECK(names.contains("shell"));
    }

    SECTION("a renamed body can still be referenced by its new name") {
        auto seq = parse(
            "load {\n"
            "    pdb tests/files/SASDJG5_single.pdb tests/files/SASDJG5_single.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
            "rename b2 core\n"
            "symmetry core c2\n"
        );
        REQUIRE(seq != nullptr);
        auto rb = seq->_get_rigidbody();
        REQUIRE(rb != nullptr);
        CHECK(rb->molecule.get_body(0).size_symmetry() == 0);
        CHECK(rb->molecule.get_body(1).size_symmetry() == 1);
    }

    SECTION("the underlying body is unaffected by the rename") {
        auto seq = parse(
            "load {\n"
            "    pdb tests/files/SASDJG5_single.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
            "rename b1 core\n"
        );
        REQUIRE(seq != nullptr);
        auto rb = seq->_get_rigidbody();
        REQUIRE(rb != nullptr);
        CHECK(rb->molecule.size_body() == 1);
        CHECK(rb->molecule.get_body(0).size_atom() > 0);
    }

    SECTION("renaming an unknown body name is rejected") {
        CHECK_THROWS(parse(
            "load {\n"
            "    pdb tests/files/SASDJG5_single.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
            "rename doesnotexist core\n"
        ));
    }

    SECTION("renaming to an already-used name is rejected") {
        CHECK_THROWS(parse(
            "load {\n"
            "    pdb tests/files/SASDJG5_single.pdb tests/files/SASDJG5_single.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
            "rename b1 b2\n"
        ));
    }

    SECTION("renaming a body to itself is rejected") {
        CHECK_THROWS(parse(
            "load {\n"
            "    pdb tests/files/SASDJG5_single.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
            "rename b1 b1\n"
        ));
    }

    SECTION("a bare rename with only one argument is rejected") {
        CHECK_THROWS(parse(
            "load {\n"
            "    pdb tests/files/SASDJG5_single.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
            "rename b1\n"
        ));
    }
}
