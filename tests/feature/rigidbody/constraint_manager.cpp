#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <rigidbody/constraints/OverlapConstraint.h>
#include <rigidbody/constraints/DistanceConstraintBond.h>
#include <rigidbody/constraints/DistanceConstraintAtom.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/Rigidbody.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <settings/All.h>

using namespace ausaxs;
using namespace ausaxs::data;
using namespace ausaxs::rigidbody;

struct fixture {
    fixture() {
        settings::molecule::implicit_hydrogens = false;
        settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::None;
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
    SECTION("Protein*") {
        constraints::ConstraintManager cm(&protein);
        CHECK(cm.molecule == &protein.molecule);
    }
}

TEST_CASE_METHOD(fixture, "ConstraintManager::add_constraint") {
    settings::general::verbose = false;
    Rigidbody protein(Molecule{ap});

    SECTION("OverlapConstraint") {
        constraints::ConstraintManager cm(&protein);
        auto initial_non_disc = cm.non_discoverable_constraints.size();
        cm.add_constraint(std::make_unique<constraints::OverlapConstraint>(&protein.molecule));
        CHECK(cm.non_discoverable_constraints.size() == initial_non_disc + 1);
    }

    SECTION("DistanceConstraintBond") {
        constraints::ConstraintManager cm(&protein);
        cm.add_constraint(std::make_unique<constraints::DistanceConstraintBond>(&protein.molecule, 0, 1));
        CHECK(cm.discoverable_constraints.size() == 1);
    }

    SECTION("Multiple") {
        constraints::ConstraintManager cm(&protein);
        auto initial_non_disc = cm.non_discoverable_constraints.size();
        cm.add_constraint(std::make_unique<constraints::OverlapConstraint>(&protein.molecule));
        cm.add_constraint(std::make_unique<constraints::DistanceConstraintBond>(&protein.molecule, 0, 1));
        cm.add_constraint(std::make_unique<constraints::DistanceConstraintBond>(&protein.molecule, 0, 2));
        CHECK(cm.non_discoverable_constraints.size() == initial_non_disc + 1);
        REQUIRE(cm.discoverable_constraints.size() == 2);
    }
}

TEST_CASE_METHOD(fixture, "ConstraintManager::evaluate") {
    settings::general::verbose = false;
    Rigidbody protein(Molecule{ap});
    SECTION("returns chi2 contribution of all constraints") {
        constraints::ConstraintManager cm(&protein);
        // Remove the automatically added OverlapConstraint so we can test distance constraints in isolation.
        cm.non_discoverable_constraints.clear();
        REQUIRE(cm.evaluate() == 0);

        cm.add_constraint(std::make_unique<constraints::DistanceConstraintBond>(&protein.molecule, 0, 1));
        auto dc1 = cm.discoverable_constraints.back().get();

        // Move body 0 toward body 1 to compress the bond (triggers non-zero evaluate)
        protein.molecule.get_body(0).translate(Vector3<double>(0, 0, 1));
        auto val = dc1->evaluate();
        CHECK(val != 0);
        CHECK(cm.evaluate() == val);

        // shift back to original position, should be zero again
        protein.molecule.get_body(0).translate(Vector3<double>(0, 0, -1));
        CHECK(dc1->evaluate() == 0);

        // shift again
        protein.molecule.get_body(0).translate(Vector3<double>(0, 0, -1));
        val = dc1->evaluate();
        CHECK(val == 0);
        CHECK(cm.evaluate() == val);

        cm.add_constraint(std::make_unique<constraints::DistanceConstraintBond>(&protein.molecule, 0, 2));
        auto dc2 = cm.discoverable_constraints.back().get();
        protein.molecule.get_body(0).translate(Vector3<double>(1, 0, 0));
        auto val2 = dc2->evaluate();
        auto val_after = dc1->evaluate();
        CHECK(val2 != 0);
        CHECK(cm.evaluate() == val_after + val2);
    }
}