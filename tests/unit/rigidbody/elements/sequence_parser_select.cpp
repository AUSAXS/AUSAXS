// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <catch2/catch_test_macros.hpp>

#include <rigidbody/sequencer/detail/SequenceParser.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/Rigidbody.h>
#include <rigidbody/selection/ManualSelect.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <settings/All.h>
#include <io/ExistingFile.h>

#include <fstream>

using namespace ausaxs;
using namespace ausaxs::rigidbody;
using namespace ausaxs::rigidbody::sequencer;

struct SequenceParserSelectFixture {
    SequenceParserSelectFixture() {
        settings::general::verbose = false;
        settings::molecule::implicit_hydrogens = false;
        settings::grid::min_bins = 250;
        settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::None;
    }

    std::unique_ptr<Sequencer> parse(const std::string& content) {
        static int counter = 0;
        std::string path = "/tmp/ausaxs_seq_select_test_" + std::to_string(counter++) + ".conf";
        std::ofstream f(path);
        f << content;
        f.close();
        SequenceParser parser;
        return parser.parse_file(path);
    }
};

TEST_CASE_METHOD(SequenceParserSelectFixture, "SequenceParser::BodySelectElement inline") {
    SECTION("selecting a body name switches to ManualSelect targeting that body") {
        auto seq = parse(
            "load {\n"
            "    pdb tests/files/SASDJG5_single.pdb tests/files/SASDJG5_single.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
            "select b2\n"
        );
        REQUIRE(seq != nullptr);
        seq->execute();
        auto rb = seq->_get_rigidbody();
        REQUIRE(rb != nullptr);
        REQUIRE(rb->body_selector != nullptr);
        CHECK(dynamic_cast<selection::ManualSelect*>(rb->body_selector.get()) != nullptr);
    }

    SECTION("selecting a known strategy name inline works without braces") {
        auto seq = parse(
            "load {\n"
            "    pdb tests/files/SASDJG5_single.pdb tests/files/SASDJG5_single.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
            "select sequential_body\n"
        );
        REQUIRE(seq != nullptr);
        seq->execute();
        auto rb = seq->_get_rigidbody();
        REQUIRE(rb != nullptr);
        REQUIRE(rb->body_selector != nullptr);
        CHECK(dynamic_cast<selection::ManualSelect*>(rb->body_selector.get()) == nullptr);
    }

    SECTION("selecting an unknown body name or strategy is rejected") {
        CHECK_THROWS(parse(
            "load {\n"
            "    pdb tests/files/SASDJG5_single.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
            "select doesnotexist\n"
        ));
    }

    SECTION("more than one inline argument is rejected") {
        CHECK_THROWS(parse(
            "load {\n"
            "    pdb tests/files/SASDJG5_single.pdb tests/files/SASDJG5_single.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
            "select b1 b2\n"
        ));
    }

    SECTION("named argument block form still works") {
        auto seq = parse(
            "load {\n"
            "    pdb tests/files/SASDJG5_single.pdb tests/files/SASDJG5_single.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
            "select {\n"
            "    point random_body\n"
            "}\n"
        );
        REQUIRE(seq != nullptr);
        seq->execute();
        auto rb = seq->_get_rigidbody();
        REQUIRE(rb != nullptr);
        REQUIRE(rb->body_selector != nullptr);
    }
}
