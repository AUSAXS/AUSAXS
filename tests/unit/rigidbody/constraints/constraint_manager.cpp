#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <data/Body.h>
#include <data/Molecule.h>
#include <rigidbody/Rigidbody.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/constraints/DistanceConstraintAtom.h>
#include <rigidbody/constraints/DistanceConstraintBond.h>
#include <rigidbody/constraints/OverlapConstraint.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/sequencer/detail/SequenceParser.h>
#include <settings/All.h>

#include <support/rb_metadata.h>
#include <support/temp_file.h>

#include <vector>

using namespace ausaxs;
using namespace ausaxs::data;
using namespace ausaxs::rigidbody;
using namespace ausaxs::rigidbody::sequencer;

struct fixture {
    fixture() {
        settings::molecule::implicit_hydrogens = false;
    }

    AtomFF a1 = AtomFF({-1, -1, -1}, form_factor::form_factor_t::C);
    AtomFF a2 = AtomFF({-1,  1, -1}, form_factor::form_factor_t::C);
    AtomFF a3 = AtomFF({-1, -1,  1}, form_factor::form_factor_t::C);
    AtomFF a4 = AtomFF({-1,  1,  1}, form_factor::form_factor_t::C);
    AtomFF a5 = AtomFF({ 1, -1, -1}, form_factor::form_factor_t::C);
    AtomFF a6 = AtomFF({ 1,  1, -1}, form_factor::form_factor_t::C);
    AtomFF a7 = AtomFF({ 1, -1,  1}, form_factor::form_factor_t::C);
    AtomFF a8 = AtomFF({ 1,  1,  1}, form_factor::form_factor_t::NH);

    Body b1 = Body(std::vector{a1, a2});
    Body b2 = Body(std::vector{a3, a4});
    Body b3 = Body(std::vector{a5, a6});
    Body b4 = Body(std::vector{a7, a8});
    std::vector<Body> ap = {b1, b2, b3, b4};
};

TEST_CASE_METHOD(fixture, "ConstraintManager::ConstraintManager") {
    settings::general::verbose = false;
    Rigidbody protein(Molecule{ap});
    test::mark_backbone_carbons(protein.molecule);
    SECTION("Protein*") {
        constraints::ConstraintManager cm(&protein);
        CHECK(cm.molecule == &protein.molecule);
    }
}

TEST_CASE_METHOD(fixture, "ConstraintManager::add_constraint") {
    settings::general::verbose = false;
    Rigidbody protein(Molecule{ap});
    test::mark_backbone_carbons(protein.molecule);

    SECTION("OverlapConstraint is rejected") {
        // the manager creates the one and only OverlapConstraint itself; a second must not be accepted
        constraints::ConstraintManager cm(&protein);
        auto initial_non_disc = cm.non_discoverable_constraints.size();
        CHECK_THROWS(cm.add_constraint(std::make_unique<constraints::OverlapConstraint>(&protein.molecule)));
        CHECK(cm.non_discoverable_constraints.size() == initial_non_disc);
    }

    SECTION("DistanceConstraintBond") {
        constraints::ConstraintManager cm(&protein);
        cm.add_constraint(std::make_unique<constraints::DistanceConstraintBond>(&protein.molecule, 0, 1));
        CHECK(cm.discoverable_constraints.size() == 1);
    }

    SECTION("DistanceConstraintAtom") {
        constraints::ConstraintManager cm(&protein);
        auto initial_non_disc = cm.non_discoverable_constraints.size();
        cm.add_constraint(std::make_unique<constraints::DistanceConstraintAtom>(&protein.molecule, 0, 0, 1, 0));
        CHECK(cm.non_discoverable_constraints.size() == initial_non_disc + 1);
    }

    SECTION("Multiple") {
        constraints::ConstraintManager cm(&protein);
        auto initial_non_disc = cm.non_discoverable_constraints.size();
        cm.add_constraint(std::make_unique<constraints::DistanceConstraintAtom>(&protein.molecule, 0, 0, 1, 0));
        cm.add_constraint(std::make_unique<constraints::DistanceConstraintBond>(&protein.molecule, 0, 1));
        cm.add_constraint(std::make_unique<constraints::DistanceConstraintBond>(&protein.molecule, 0, 1));
        CHECK(cm.non_discoverable_constraints.size() == initial_non_disc + 1);
        REQUIRE(cm.discoverable_constraints.size() == 2);
    }
}

TEST_CASE_METHOD(fixture, "ConstraintManager::evaluate") {
    settings::general::verbose = false;
    Rigidbody protein(Molecule{ap});
    test::mark_backbone_carbons(protein.molecule);
    SECTION("returns chi2 contribution of all constraints") {
        constraints::ConstraintManager cm(&protein);
        // Remove the automatically added OverlapConstraint so we can test distance constraints in isolation.
        cm.non_discoverable_constraints.clear();
        REQUIRE(cm.evaluate() == 0);

        cm.add_constraint(std::make_unique<constraints::DistanceConstraintBond>(&protein.molecule, 0, 1));
        auto* dc1 = cm.discoverable_constraints.back().get();

        // Move body 0 toward body 1 to compress the bond (triggers non-zero evaluate)
        protein.molecule.get_body(0).translate(Vector3<double>(0, 0, 1));
        auto val = dc1->evaluate();
        CHECK(val != 0);
        CHECK(cm.evaluate() == val);

        // shift back to original position, should be zero again
        protein.molecule.get_body(0).translate(Vector3<double>(0, 0, -1));
        CHECK(dc1->evaluate() == 0);

        // shift again (body is now further than d_target, so still penalized)
        protein.molecule.get_body(0).translate(Vector3<double>(0, 0, -1));
        val = dc1->evaluate();
        CHECK(val != 0);
        CHECK(cm.evaluate() == val);

        cm.add_constraint(std::make_unique<constraints::DistanceConstraintBond>(&protein.molecule, 0, 1));
        auto* dc2 = cm.discoverable_constraints.back().get();
        protein.molecule.get_body(0).translate(Vector3<double>(1, 0, 0));
        auto val2 = dc2->evaluate();
        auto val_after = dc1->evaluate();
        CHECK(val2 != 0);
        CHECK(cm.evaluate() == val_after + val2);
    }
}

namespace {
    void common_settings() {
        settings::general::verbose = false;
        settings::molecule::implicit_hydrogens = false;
        settings::grid::min_bins = 25; // let the grid size itself to the conformation the optimization wanders into
    }
}

// The per-body constraint map is keyed by body index, so a body added after the manager was built has no entry of its own. Every selector reads the map each
// iteration, so a stale one aborts the optimization outright rather than merely answering wrongly.
TEST_CASE("ConstraintManager: the per-body map covers bodies added after construction") {
    common_settings();

    std::vector<Body> bodies;
    bodies.reserve(2);
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
TEST_CASE("ConstraintManager: an unconstrained split optimizes without tripping the per-body map") {
    common_settings();

    test::TempFile config(".conf",
        "load {\n"
        "    pdb tests/files/LAR1-4.pdb\n"
        "    saxs tests/files/LAR1-2.dat\n"
        "}\n"
        "split b1 99 202\n"
        "parameters {\n"
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
