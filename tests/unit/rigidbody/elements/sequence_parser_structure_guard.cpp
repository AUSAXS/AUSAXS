// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <rigidbody/sequencer/detail/SequenceParser.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/Rigidbody.h>
#include <data/Molecule.h>
#include <settings/All.h>
#include <io/ExistingFile.h>

#include <support/temp_file.h>

#include <string>

using namespace ausaxs;
using namespace ausaxs::rigidbody;
using namespace ausaxs::rigidbody::sequencer;

namespace {
    std::unique_ptr<Sequencer> parse(const std::string& content) {
        test::TempFile config("ausaxs_seq_structure_guard_test", ".conf", content);
        SequenceParser parser;
        return parser.parse_file(config);
    }

    // One chain cut into five backbone-adjacent bodies, so "autoconstrain backbone" actually produces constraints to guard against.
    std::string load_split() {
        return
            "load {\n"
            "    pdb tests/files/LAR1-4.pdb\n"
            "    saxs tests/files/LAR1-2.dat\n"
            "    split 9 99 202 292\n"
            "}\n";
    }
}

// Constraints and the symmetry target pool are both indexed by body index, and constraints are never reindexed. Changing the set of bodies afterwards would
// silently leave both pointing at the wrong bodies, so it is rejected outright.
TEST_CASE("SequenceParser: the body set cannot be changed once constraints exist") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::grid::min_bins = 25; // let the grid size itself: "copy" places its clone 2*Rg away, which overflows a fixed small grid
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::None;

    // the guard runs before the element does any work of its own, so every mutating element can be checked the same way
    auto mutation = GENERATE(as<std::string>{},
        "delete b2\n",
        "merge b1 b2\n",
        "copy b6 b1\n",
        "split b1 5\n",
        "convert_to_symmetry {\n    type c4\n    bodies b1 b2 b3 b4\n    tolerance 100\n}\n"
    );

    // sanity: the fixture really does declare constraints, otherwise the guard would have nothing to trip on and the test would pass vacuously
    auto seq = parse(load_split() + "autoconstrain backbone\n");
    REQUIRE_FALSE(seq->_get_rigidbody()->constraints->discoverable_constraints.empty());

    CHECK_THROWS(parse(load_split() + "autoconstrain backbone\n" + mutation));
}

TEST_CASE("SequenceParser: the body set may be changed while no constraints exist") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::grid::min_bins = 25; // let the grid size itself: "copy" places its clone 2*Rg away, which overflows a fixed small grid
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::None;

    auto mutation = GENERATE(as<std::string>{},
        "delete b2\n",
        "merge b1 b2\n",
        "copy b6 b1\n",
        "split b1 5\n"
    );

    SECTION("with no constraints at all") {
        CHECK_NOTHROW(parse(load_split() + mutation));
    }

    SECTION("declared afterwards, which is the supported ordering") {
        CHECK_NOTHROW(parse(load_split() + mutation + "autoconstrain backbone\n"));
    }
}
