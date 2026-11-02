// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <io/ExistingFile.h>
#include <rigidbody/Rigidbody.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/sequencer/detail/SequenceParser.h>
#include <rigidbody/sequencer/detail/parse_error.h>
#include <settings/All.h>

#include <support/temp_file.h>

#include <string>

using namespace ausaxs;
using namespace ausaxs::rigidbody;
using namespace ausaxs::rigidbody::sequencer;

namespace {
    std::unique_ptr<Sequencer> parse(const std::string& content) {
        test::TempFile config(".conf", content);
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

// The guard above is one rung of a general ordering: load, then changes to the set of bodies, then constraints, then the optimization.
// Elements which only need a loaded structure may appear anywhere after it, and elements which use nothing at parse time anywhere at all.
TEST_CASE("SequenceParser: elements must appear in phase order") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::grid::min_bins = 25;

    auto matches = [] (const std::string& script, const std::string& expected) {
        try {
            parse(script);
        } catch (const sequencer::except::parse_error& e) {
            return std::string(e.what()).find(expected) != std::string::npos;
        }
        return false;
    };

    SECTION("an element using the structure cannot come before the load") {
        // each of these reaches for the structure at parse time, so it must be rejected before it gets that far
        auto element = GENERATE(as<std::string>{},
            "parameters {\n    iterations 10\n    translate 1\n}\n",
            "select random_body\n",
            "transform rigid\n",
            "rename b1 first\n",
            "relative_hydration b1 0.5\n",
            "symmetry c2\n",
            "delete b2\n",
            "split b1 5\n",
            "copy b6 b1\n",
            "merge b1 b2\n",
            "constrain {\n    body1 b1\n    body2 b2\n    type bond\n}\n",
            "autoconstrain backbone\n",
            "loop 1\nend\n",
            "optimize_once\nend\n"
        );
        CHECK(matches(element + load_split(), "at line 1: It requires a loaded structure"));
    }

    SECTION("elements which only act when they run may come before the load") {
        CHECK_NOTHROW(parse("seed 1\nprint \"hello\"\nsave out.pdb\noutput somewhere/\n" + load_split()));
    }

    SECTION("a script must load a structure") {
        CHECK_THROWS_AS(parse("seed 1\nprint \"hello\"\n"), sequencer::except::parse_error);
    }

    SECTION("only a single load is allowed") {
        CHECK(matches(load_split() + load_split(), "at line 6: Only one \"load\" element is allowed per script, and the structure was already loaded by \"load\" at line 1."));
    }

    SECTION("the body set cannot change after constraints, and the error names the culprit") {
        CHECK(matches(load_split() + "autoconstrain backbone\ndelete b2\n", "\"delete\" at line 7: It must come before \"autoconstrain\" at line 6."));
    }

    SECTION("nothing can be set up once the optimization has started") {
        auto element = GENERATE(as<std::string>{},
            "delete b2\n",
            "symmetry c2\n",
            "autoconstrain backbone\n",
            "constrain {\n    body1 b1\n    body2 b2\n    type bond\n}\n"
        );
        CHECK(matches(load_split() + "loop 1\nend\n" + element, "It must come before \"loop\" at line 6, which starts the optimization."));

        // the setup is closed by the first optimization element even when it is nested inside a loop
        CHECK(matches(load_split() + "loop 1\n    optimize_once\n    end\n    " + element + "end\n", "It must come before \"loop\" at line 6"));
    }

    SECTION("elements which only need the structure may follow the constraints and the optimization") {
        CHECK_NOTHROW(parse(load_split() +
            "autoconstrain backbone\n"
            "parameters {\n    iterations 10\n    translate 1\n}\n"
            "select random_body\n"
            "loop 1\n"
            "    transform single\n"
            "    optimize_once\n"
            "    end\n"
            "end\n"
            "rename b1 first\n"
            "parameters {\n    iterations 10\n    rotate 1\n}\n"
            "loop\nend\n"
        ));
    }
}
