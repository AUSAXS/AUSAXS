// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <catch2/catch_test_macros.hpp>

#include <rigidbody/sequencer/detail/SequenceParser.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/Rigidbody.h>
#include <fitter/FitResult.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <settings/All.h>

#include <support/temp_file.h>

#include <vector>

using namespace ausaxs;
using namespace ausaxs::data;
using namespace ausaxs::rigidbody;
using namespace ausaxs::rigidbody::sequencer;

namespace {
    void common_settings() {
        settings::general::verbose = false;
        settings::molecule::implicit_hydrogens = false;
        settings::grid::min_bins = 25; // let the grid size itself to the conformation the optimization wanders into
        settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::None;
    }
}

// The per-body constraint map is keyed by body index, so a body added after the manager was built has no entry of its own. Every selector reads the map each
// iteration, so a stale one aborts the optimization outright rather than merely answering wrongly.
TEST_CASE("ConstraintManager: the per-body map covers bodies added after construction") {
    common_settings();

    std::vector<Body> bodies;
    for (int i = 0; i < 2; ++i) {bodies.emplace_back(std::vector{AtomFF({5.0*i, 0, 0}, form_factor::form_factor_t::C)});}
    Rigidbody rb{Molecule{std::move(bodies)}};
    REQUIRE(rb.constraints->discoverable_constraints.empty());

    SECTION("an unconstrained body still has an entry, and it is empty") {
        CHECK(rb.constraints->get_body_constraints(0).empty());
        CHECK(rb.constraints->get_body_constraints(1).empty());
    }

    SECTION("a body gained after the map was first built is covered once it is invalidated") {
        REQUIRE(rb.constraints->get_body_constraints(0).empty()); // force the map to be built while the molecule still has two bodies
        rb.molecule.get_bodies().emplace_back(std::vector{AtomFF({10, 0, 0}, form_factor::form_factor_t::C)});
        rb.constraints->invalidate();
        CHECK(rb.constraints->get_body_constraints(2).empty());
    }
}

// Regression: splitting a body hands the optimizer indices the constraint map predates. Splitting before declaring constraints is the supported ordering - the
// structure guard forbids only the reverse - so the map has to catch up on its own instead of the selector walking off the end of it.
TEST_CASE("ConstraintManager: an unconstrained split optimizes without tripping the per-body map", "[files]") {
    common_settings();

    test::TempFile config(".conf",
        "load {\n"
        "    pdb tests/files/LAR1-4.pdb\n"
        "    saxs tests/files/LAR1-2.dat\n"
        "}\n"
        "split b1 99 202\n"
        "parameter_strategy {\n"
        "    iterations 3\n"
        "    translate 1\n"
        "    rotate 1\n"
        "}\n"
        "loop 3\n"
        "    optimize_once\n"
        "    end\n"
        "end\n"
    );

    SequenceParser parser;
    auto seq = parser.parse_file(config);
    REQUIRE(seq != nullptr);

    // sanity: the scenario really is the unconstrained one, and the split really did grow the body set
    REQUIRE(seq->_get_rigidbody()->constraints->discoverable_constraints.empty());
    REQUIRE(seq->_get_molecule()->size_body() == 3);

    CHECK_NOTHROW(seq->execute());
}
