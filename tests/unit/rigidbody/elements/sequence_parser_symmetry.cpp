// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <catch2/catch_test_macros.hpp>

#include <rigidbody/sequencer/detail/SequenceParser.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/Rigidbody.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/constraints/IDistanceConstraint.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <settings/All.h>
#include <io/ExistingFile.h>

#include <fstream>

using namespace ausaxs;
using namespace ausaxs::rigidbody;
using namespace ausaxs::rigidbody::sequencer;
using namespace ausaxs::rigidbody::constraints;

struct SequenceParserSymmetryFixture {
    SequenceParserSymmetryFixture() {
        settings::general::verbose = false;
        settings::molecule::implicit_hydrogens = false;
        settings::grid::min_bins = 250;
        settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::None;
    }

    std::unique_ptr<Sequencer> parse(const std::string& content) {
        static int counter = 0;
        std::string path = "/tmp/ausaxs_seq_sym_test_" + std::to_string(counter++) + ".conf";
        std::ofstream f(path);
        f << content;
        f.close();
        SequenceParser parser;
        return parser.parse(io::ExistingFile(path));
    }
};

TEST_CASE_METHOD(SequenceParserSymmetryFixture, "SequenceParser::SymmetryElement") {
    SECTION("c2 symmetry creates one symmetry on the body") {
        auto seq = parse(
            "load {\n"
            "    pdb tests/files/SASDJG5_single.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
            "symmetry c2\n"
        );
        REQUIRE(seq != nullptr);
        auto rb = seq->_get_rigidbody();
        REQUIRE(rb != nullptr);
        CHECK(rb->molecule.get_body(0).size_symmetry() == 1);
    }

    SECTION("no symmetry directive leaves body without symmetries") {
        auto seq = parse(
            "load {\n"
            "    pdb tests/files/SASDJG5_single.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
        );
        REQUIRE(seq != nullptr);
        auto rb = seq->_get_rigidbody();
        REQUIRE(rb != nullptr);
        CHECK(rb->molecule.get_body(0).size_symmetry() == 0);
    }

    SECTION("two symmetry directives produce two symmetries on the body") {
        auto seq = parse(
            "load {\n"
            "    pdb tests/files/SASDJG5_single.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
            "symmetry c2\n"
            "symmetry c2\n"
        );
        REQUIRE(seq != nullptr);
        auto rb = seq->_get_rigidbody();
        REQUIRE(rb != nullptr);
        CHECK(rb->molecule.get_body(0).size_symmetry() == 2);
    }

    SECTION("symmetry applied to one body does not affect other bodies") {
        auto seq = parse(
            "load {\n"
            "    pdb tests/files/SASDJG5_single.pdb tests/files/SASDJG5_single.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
            "symmetry {\n"
            "    b2 c2\n"
            "}\n"
        );
        REQUIRE(seq != nullptr);
        auto rb = seq->_get_rigidbody();
        REQUIRE(rb != nullptr);
        CHECK(rb->molecule.get_body(0).size_symmetry() == 0);
        CHECK(rb->molecule.get_body(1).size_symmetry() == 1);
    }
}

TEST_CASE_METHOD(SequenceParserSymmetryFixture, "SequenceParser::ConstraintElement real-real") {
    auto seq = parse(
        "load {\n"
        "    pdb tests/files/SASDJG5_single.pdb tests/files/SASDJG5_single.pdb\n"
        "    saxs tests/files/SASDJG5.dat\n"
        "}\n"
        "constrain {\n"
        "    body1 b1\n"
        "    body2 b2\n"
        "    type attract\n"
        "    distance 30\n"
        "}\n"
    );
    REQUIRE(seq != nullptr);
    auto rb = seq->_get_rigidbody();
    REQUIRE(rb != nullptr);
    // non_discoverable_constraints[0] is the pre-added OverlapConstraint; ours is at the back
    REQUIRE(rb->constraints->non_discoverable_constraints.size() >= 2);
    auto* c = dynamic_cast<IDistanceConstraint*>(rb->constraints->non_discoverable_constraints.back().get());
    REQUIRE(c != nullptr);

    SECTION("ibody1 and ibody2 reference the two distinct bodies") {
        CHECK(c->ibody1 == 0);
        CHECK(c->ibody2 == 1);
    }

    SECTION("both isym values indicate real (non-symmetry) bodies") {
        CHECK(c->isym1 == std::make_pair(-1, 0));
        CHECK(c->isym2 == std::make_pair(-1, 0));
    }
}

TEST_CASE_METHOD(SequenceParserSymmetryFixture, "SequenceParser::ConstraintElement real-symmetry") {
    auto seq = parse(
        "load {\n"
        "    pdb tests/files/SASDJG5_single.pdb\n"
        "    saxs tests/files/SASDJG5.dat\n"
        "}\n"
        "symmetry c2\n"
        "constrain {\n"
        "    body1 b1s1\n"
        "    body2 b1\n"
        "    type attract\n"
        "    distance 30\n"
        "}\n"
    );
    REQUIRE(seq != nullptr);
    auto rb = seq->_get_rigidbody();
    REQUIRE(rb != nullptr);
    REQUIRE(rb->constraints->non_discoverable_constraints.size() >= 2);
    auto* c = dynamic_cast<IDistanceConstraint*>(rb->constraints->non_discoverable_constraints.back().get());
    REQUIRE(c != nullptr);

    SECTION("both constrained atoms belong to body 0") {
        CHECK(c->ibody1 == 0);
        CHECK(c->ibody2 == 0);
    }

    SECTION("body1 isym tracks the first symmetry, first replica") {
        CHECK(c->isym1 == std::make_pair(0, 1));
    }

    SECTION("body2 isym indicates a real (non-symmetry) body") {
        CHECK(c->isym2 == std::make_pair(-1, 0));
    }
}

TEST_CASE_METHOD(SequenceParserSymmetryFixture, "SequenceParser::ConstraintElement symmetry-symmetry") {
    auto seq = parse(
        "load {\n"
        "    pdb tests/files/SASDJG5_single.pdb tests/files/SASDJG5_single.pdb\n"
        "    saxs tests/files/SASDJG5.dat\n"
        "}\n"
        "symmetry {\n"
        "    b1 c2\n"
        "    b2 c2\n"
        "}\n"
        "constrain {\n"
        "    body1 b1s1\n"
        "    body2 b2s1\n"
        "    type attract\n"
        "    distance 50\n"
        "}\n"
    );
    REQUIRE(seq != nullptr);
    auto rb = seq->_get_rigidbody();
    REQUIRE(rb != nullptr);
    REQUIRE(rb->constraints->non_discoverable_constraints.size() >= 2);
    auto* c = dynamic_cast<IDistanceConstraint*>(rb->constraints->non_discoverable_constraints.back().get());
    REQUIRE(c != nullptr);

    SECTION("body indices reference the two distinct bodies") {
        CHECK(c->ibody1 == 0);
        CHECK(c->ibody2 == 1);
    }

    SECTION("both isym values track the first symmetry, first replica of their respective bodies") {
        CHECK(c->isym1 == std::make_pair(0, 1));
        CHECK(c->isym2 == std::make_pair(0, 1));
    }
}
