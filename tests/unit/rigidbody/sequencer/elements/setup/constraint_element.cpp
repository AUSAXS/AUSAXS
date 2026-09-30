#include <catch2/catch_test_macros.hpp>

#include <data/Body.h>
#include <data/Molecule.h>
#include <data/symmetry/CompositeSymmetry.h>
#include <data/symmetry/ReferenceSymmetry.h>
#include <fitter/FitResult.h>  // IWYU pragma: keep
#include <io/ExistingFile.h>
#include <io/Folder.h>
#include <math/Vector3.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/constraints/IDistanceConstraint.h>
#include <rigidbody/parameters/UniformParameterGenerator.h>
#include <rigidbody/Rigidbody.h>
#include <rigidbody/sequencer/detail/SequenceParser.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/transform/TransformStrategy.h>  // IWYU pragma: keep
#include <settings/All.h>

#include <support/temp_file.h>

#include <algorithm>

using namespace ausaxs;
using namespace ausaxs::data;
using namespace ausaxs::rigidbody;
using namespace ausaxs::rigidbody::sequencer;
using namespace ausaxs::rigidbody::constraints;

struct SequencerElementsFixture {
    SequencerElementsFixture() {
        settings::general::verbose = false;
        settings::molecule::implicit_hydrogens = false;
        settings::grid::min_bins = 250;
    }
};

namespace {
    std::unique_ptr<Sequencer> parse_sequence(const std::string& script) {
        SequenceParser parser;
        return parser.parse_text(script + (script.find("\nloop ") != std::string::npos ? "end\n" : ""));
    }
}

TEST_CASE_METHOD(SequencerElementsFixture, "SequencerElements::ConstraintElement") {
    SECTION("Add distance constraint center mass") {
        REQUIRE_NOTHROW(
            parse_sequence(
                "load {\n"
                "    pdb tests/files/LAR1-2.pdb\n"
                "    saxs tests/files/LAR1-2.dat\n"
                "    split 9 99\n"
                "}\n"
                "constrain {\n"
                "    first b1\n"
                "    second b2\n"
                "    type cm\n"
                "}\n"
            )
        );
    }

    SECTION("Add distance constraint closest") {
        REQUIRE_NOTHROW(
            parse_sequence(
                "load {\n"
                "    pdb tests/files/LAR1-2.pdb\n"
                "    saxs tests/files/LAR1-2.dat\n"
                "    split 9 99\n"
                "}\n"
                "constrain {\n"
                "    first b1\n"
                "    second b2\n"
                "    type bond\n"
                "}\n"
            )
        );
    }
}

// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

struct SequenceParserSymmetryFixture {
    SequenceParserSymmetryFixture() {
        settings::general::verbose = false;
        settings::molecule::implicit_hydrogens = false;
        settings::grid::min_bins = 250;
    }

    static std::unique_ptr<Sequencer> parse(const std::string& content) {
        test::TempFile config(".conf", content);
        SequenceParser parser;
        return parser.parse_file(config);
    }
};

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
    auto* rb = seq->_get_rigidbody();
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
    auto* rb = seq->_get_rigidbody();
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
        "symmetry b1 c2\n"
        "symmetry b2 c2\n"
        "constrain {\n"
        "    body1 b1s1\n"
        "    body2 b2s1\n"
        "    type attract\n"
        "    distance 50\n"
        "}\n"
    );
    REQUIRE(seq != nullptr);
    auto* rb = seq->_get_rigidbody();
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
