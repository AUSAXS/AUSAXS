#include <catch2/catch_test_macros.hpp>

#include <rigidbody/constraints/OverlapConstraint.h>
#include <rigidbody/constraints/DistanceConstraintBond.h>
#include <rigidbody/constraints/DistanceConstraintAtom.h>
#include <rigidbody/constraints/DistanceConstraintCM.h>
#include <rigidbody/constraints/AttractorConstraint.h>
#include <rigidbody/constraints/RepellerConstraint.h>
#include <rigidbody/constraints/FixedConstraint.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <settings/All.h>

using namespace ausaxs;
using namespace ausaxs::data;
using namespace ausaxs::rigidbody::constraints;

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

TEST_CASE_METHOD(fixture, "Constraints::basic_evaluate") {
    Molecule mol(ap);

    SECTION("DistanceConstraintBond") {
        DistanceConstraintBond c(&mol, 0, 1);
        CHECK(c.evaluate() == 0);
        // move body0 toward body1
        mol.get_body(0).translate(Vector3<double>(0,0,1));
        CHECK(c.evaluate() > 0);
    }

    SECTION("DistanceConstraintAtom") {
        // choose atom indices guaranteed to exist (first carbon atoms)
        DistanceConstraintAtom c(&mol, 0, 0, 1, 0);
        CHECK(c.evaluate() == 0);
        mol.get_body(0).translate(Vector3<double>(0,0,1));
        CHECK(c.evaluate() > 0);
    }

    SECTION("DistanceConstraintCM") {
        DistanceConstraintCM c(&mol, 0, 1);
        CHECK(c.evaluate() == 0);
        mol.get_body(0).translate(Vector3<double>(0,0,1));
        CHECK(c.evaluate() > 0);
    }

    SECTION("AttractorConstraint and RepellerConstraint") {
        // compute current CM-atom distance used by DistanceConstraintCM
        auto atom1 = mol.get_body(0).get_atom(0).coordinates();
        auto atom2 = mol.get_body(1).get_atom(0).coordinates();
        double current = atom1.distance(atom2);

        AttractorConstraint a(&mol, current - 0.5, 0, 1);
        // since target < current, attractor should penalize (distance > target)
        CHECK(a.evaluate() > 0);

        RepellerConstraint r(&mol, current + 0.5, 0, 1);
        // since target > current, repeller should penalize (distance < target)
        CHECK(r.evaluate() > 0);
    }

    SECTION("OverlapConstraint") {
        OverlapConstraint o(&mol);
        // initialize() runs in ctor; evaluate should be non-negative
        double v0 = o.evaluate();
        CHECK(v0 >= 0);
        mol.get_body(0).translate(Vector3<double>(0,0,1));
        double v1 = o.evaluate();
        CHECK(v1 >= 0);
        // values may change after translation
        CHECK(v1 != v0);
    }

    SECTION("FixedConstraint") {
        FixedConstraint f(&mol, 0, 1);
        CHECK(f.evaluate() == 0);
    }
}
