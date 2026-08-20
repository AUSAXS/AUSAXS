#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <rigidbody/constraints/OverlapConstraint.h>
#include <rigidbody/constraints/DistanceConstraintBond.h>
#include <rigidbody/constraints/DistanceConstraintAtom.h>
#include <rigidbody/constraints/DistanceConstraintCM.h>
#include <rigidbody/constraints/AttractorConstraint.h>
#include <rigidbody/constraints/RepellerConstraint.h>
#include <rigidbody/constraints/FixedConstraint.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <data/symmetry/PointSymmetry.h>
#include <math/MatrixUtils.h>
#include <settings/All.h>

#include <support/rb_metadata.h>

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
    test::mark_backbone_carbons(mol);

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

// A constraint may target one of a body's symmetry copies rather than the body itself. The copies are anchored on the body's centre of mass, so the
// constraint has to place them from that centre - not from the representative atom it happens to hold - and it has to keep doing so once the body has
// been rotated away from the frame its symmetry was defined in.
TEST_CASE("DistanceConstraint: a symmetry replica is placed on the body's centre of mass") {
    settings::molecule::implicit_hydrogens = false;
    settings::molecule::center = false;

    // exposes the protected distance evaluation so the constraint can be re-checked after the body moves
    struct Probe : DistanceConstraintAtom {
        using DistanceConstraintAtom::DistanceConstraintAtom;
        using DistanceConstraintAtom::evaluate_distance;
    };

    // deliberately off-centre and asymmetric, so anchoring on the wrong point cannot coincidentally give the right answer
    std::vector<AtomFF> body_atoms = {
        AtomFF({ 1.0,  0.0,  0.0}, form_factor::form_factor_t::C),
        AtomFF({-1.0,  0.5,  0.2}, form_factor::form_factor_t::C),
        AtomFF({ 0.3, -1.4,  0.9}, form_factor::form_factor_t::C),
        AtomFF({ 2.1,  0.7, -1.1}, form_factor::form_factor_t::C)
    };
    Molecule mol({Body{body_atoms}, Body{std::vector<AtomFF>{AtomFF({20, 0, 0}, form_factor::form_factor_t::C)}}});
    mol.get_body(0).symmetry().add(std::make_unique<symmetry::PointSymmetry>(
        Vector3<double>{8, 0, 0}, Vector3<double>{0.3, 0.7, 0.2} // a rotated copy: a wrong anchor drops the rotation entirely
    ));

    // ground truth: the materialised assembly, laid out as [body | copy 1]
    auto explicit_distance = [&mol] () {
        auto structure = mol.get_body(0).symmetry().explicit_structure();
        return structure.atoms[mol.get_body(0).size_atom()].coordinates().distance(mol.get_body(1).get_atom(0).coordinates());
    };

    Probe c(&mol, 0, 0, 1, 0, {0, 1}, {-1, -1}); // body 0, atom 0, against copy 1 of its symmetry 0
    REQUIRE_THAT(c.d_target, Catch::Matchers::WithinAbs(explicit_distance(), 1e-9));

    SECTION("and follows the body through a rigid motion") {
        auto R = matrix::rotation_matrix<double>({0.4, -0.9, 0.3});
        auto pivot = mol.get_body(0).get_cm();
        mol.get_body(0).translate(-pivot);
        mol.get_body(0).rotate(R);
        mol.get_body(0).translate(pivot + Vector3<double>{3, -2, 5});
        mol.get_body(0).symmetry().set_orientation(R);

        CHECK_THAT(c.evaluate_distance(0, 0), Catch::Matchers::WithinAbs(explicit_distance(), 1e-9));
    }
}
