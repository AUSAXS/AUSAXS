#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <data/Body.h>
#include <data/Molecule.h>
#include <data/symmetry/CyclicSymmetry.h>
#include <data/symmetry/PointSymmetry.h>
#include <data/symmetry/PredefinedSymmetries.h>
#include <math/MatrixUtils.h>
#include <rigidbody/constraints/DistanceConstraintAtom.h>
#include <settings/All.h>

using namespace ausaxs;
using namespace ausaxs::data;
using namespace ausaxs::rigidbody;

struct fixture {
    AtomFF a1 = AtomFF({-1, -1, -1}, form_factor::form_factor_t::C);
    AtomFF a2 = AtomFF({-1,  1, -1}, form_factor::form_factor_t::C);
    AtomFF a3 = AtomFF({-1, -1,  1}, form_factor::form_factor_t::C);
    AtomFF a4 = AtomFF({-1,  1,  1}, form_factor::form_factor_t::C);
    AtomFF a5 = AtomFF({ 1, -1, -1}, form_factor::form_factor_t::C);
    AtomFF a6 = AtomFF({ 1,  1, -1}, form_factor::form_factor_t::C);
    AtomFF a7 = AtomFF({ 1, -1,  1}, form_factor::form_factor_t::C);
    AtomFF a8 = AtomFF({ 1,  1,  1}, form_factor::form_factor_t::NH);

    Body b1 = Body(std::vector<AtomFF>{a1, a2});
    Body b2 = Body(std::vector<AtomFF>{a3, a4});
    Body b3 = Body(std::vector<AtomFF>{a5, a6});
    Body b4 = Body(std::vector<AtomFF>{a7, a8});
    std::vector<Body> ap = {b1, b2, b3, b4};
};

TEST_CASE_METHOD(fixture, "DistanceConstraintAtom::constructor") {
    settings::molecule::implicit_hydrogens = false;
    Molecule protein = Molecule(ap);

    SECTION("basic constructor") {
        constraints::DistanceConstraintAtom c(&protein, 0, 0, 1, 1);
        CHECK(c.ibody1 == 0);
        CHECK(c.ibody2 == 1);
        CHECK(c.iatom1 == 0);
        CHECK(c.iatom2 == 1);
        CHECK(c.isym1 == std::make_pair(-1, -1));
        CHECK(c.isym2 == std::make_pair(-1, -1));

        CHECK(c.get_body1() == b1);
        CHECK(c.get_body2() == b2);
        CHECK(c.get_atom1() == a1);
        CHECK(c.get_atom2() == a4);
    }

    SECTION("atom-reference constructor") {
        constraints::DistanceConstraintAtom c(&protein, a1, a4);
        CHECK(c.ibody1 == 0);
        CHECK(c.ibody2 == 1);
        CHECK(c.iatom1 == 0);
        CHECK(c.iatom2 == 1);

        CHECK(c.get_body1() == b1);
        CHECK(c.get_body2() == b2);
        CHECK(c.get_atom1() == a1);
        CHECK(c.get_atom2() == a4);
    }

    SECTION("constructor stores symmetry indices") {
        protein.get_body(0).symmetry().add(symmetry::type::c2);
        constraints::DistanceConstraintAtom c(&protein, 0, 0, 1, 0, {0, 1}, {-1, -1});
        CHECK(c.isym1 == std::make_pair(0, 1));
        CHECK(c.isym2 == std::make_pair(-1, -1));
    }

    SECTION("throws with non-carbon atoms") {
        CHECK_THROWS(constraints::DistanceConstraintAtom(&protein, a1, a8));
    }

    SECTION("throws within same body") {
        CHECK_THROWS(constraints::DistanceConstraintAtom(&protein, a1, a2));
    }
}

TEST_CASE_METHOD(fixture, "DistanceConstraintAtom::evaluate") {
    settings::molecule::implicit_hydrogens = false;
    Molecule protein = Molecule(ap);

    SECTION("relaxed") {
        constraints::DistanceConstraintAtom c13(&protein, 0, 0, 1, 0);
        constraints::DistanceConstraintAtom c24(&protein, 0, 1, 1, 1);
        CHECK(c13.evaluate() == 0);
        CHECK(c24.evaluate() == 0);
    }

    SECTION("stretched") {
        constraints::DistanceConstraintAtom c13(&protein, 0, 0, 1, 0);
        constraints::DistanceConstraintAtom c24(&protein, 0, 1, 1, 1);

        protein.get_body(0).translate(Vector3<double>(1, 0, 0));
        CHECK(c13.evaluate() != 0);
        CHECK(c24.evaluate() != 0);

        protein.get_body(0).translate(Vector3<double>(-1, 0, 0));
        CHECK(c13.evaluate() == 0);
        CHECK(c24.evaluate() == 0);
    }

    SECTION("only affected body changes penalty") {
        constraints::DistanceConstraintAtom c13(&protein, 0, 0, 1, 0);
        constraints::DistanceConstraintAtom c35(&protein, 1, 0, 2, 0);

        protein.get_body(0).translate(Vector3<double>(1, 0, 0));
        CHECK(c13.evaluate() != 0);
        CHECK(c35.evaluate() == 0);

        // Restore body0 and translate body1: both constraints are affected
        protein.get_body(0).translate(Vector3<double>(-1, 0, 0));
        protein.get_body(1).translate(Vector3<double>(1, 0, 0));
        CHECK(c13.evaluate() != 0);
        CHECK(c35.evaluate() != 0);
    }

    SECTION("symmetries on unrelated body do not affect real-real constraint") {
        constraints::DistanceConstraintAtom c(&protein, 0, 0, 1, 0);
        protein.get_body(2).symmetry().add(symmetry::type::c2);
        auto* sym2 = static_cast<symmetry::CyclicSymmetry*>(protein.get_body(2).symmetry().get(0));
        sym2->_initial_relation.translation = {3, 0, 0};
        CHECK(c.evaluate() == 0);

        protein.get_body(2).translate(Vector3<double>(5, 0, 0));
        CHECK(c.evaluate() == 0);
    }
}

TEST_CASE_METHOD(fixture, "DistanceConstraintAtom::evaluate with symmetry") {
    settings::molecule::implicit_hydrogens = false;
    Molecule protein = Molecule(ap);

    // c2 predefined: rotation by pi around z, zero initial translation.
    // Set a non-zero initial translation directly so the symmetric replica
    // is at a distinct position from the original atom.
    protein.get_body(0).symmetry().add(symmetry::type::c2);
    auto* sym = static_cast<symmetry::CyclicSymmetry*>(protein.get_body(0).symmetry().get(0));
    sym->_initial_relation.translation = {5, 0, 0};

    SECTION("relaxed") {
        // isym1 = {0, 1}: symmetry #0 of body0, replica 1
        constraints::DistanceConstraintAtom c(&protein, 0, 0, 1, 0, {0, 1}, {-1, -1});
        CHECK(c.evaluate() == 0);
    }

    SECTION("stretched") {
        constraints::DistanceConstraintAtom c(&protein, 0, 0, 1, 0, {0, 1}, {-1, -1});

        protein.get_body(0).translate(Vector3<double>(1, 0, 0));
        CHECK(c.evaluate() != 0);

        protein.get_body(0).translate(Vector3<double>(-1, 0, 0));
        CHECK(c.evaluate() == 0);
    }

    SECTION("unrelated body translation does not affect constraint") {
        constraints::DistanceConstraintAtom c(&protein, 0, 0, 1, 0, {0, 1}, {-1, -1});

        protein.get_body(2).translate(Vector3<double>(1, 0, 0));
        CHECK(c.evaluate() == 0);
    }

    SECTION("changing tracked symmetry translation changes result") {
        constraints::DistanceConstraintAtom c(&protein, 0, 0, 1, 0, {0, 1}, {-1, -1});
        sym->_initial_relation.translation = {7, 0, 0}; // was 5
        CHECK(c.evaluate() != 0);
        sym->_initial_relation.translation = {5, 0, 0};
        CHECK(c.evaluate() == 0);
    }

    SECTION("changing untracked symmetry on constrained body does not affect result") {
        protein.get_body(0).symmetry().add(symmetry::type::c2); // sym index 1, not tracked by isym1={0,1}
        auto* sym1 = static_cast<symmetry::CyclicSymmetry*>(protein.get_body(0).symmetry().get(1));
        sym1->_initial_relation.translation = {3, 0, 0};
        constraints::DistanceConstraintAtom c(&protein, 0, 0, 1, 0, {0, 1}, {-1, -1});
        sym1->_initial_relation.translation = {9, 0, 0}; // modify sym1 - constraint tracks sym0
        CHECK(c.evaluate() == 0);
    }
}

TEST_CASE_METHOD(fixture, "DistanceConstraintAtom::evaluate symmetry-symmetry") {
    settings::molecule::implicit_hydrogens = false;
    Molecule protein = Molecule(ap);

    // Both body0 and body1 get a c2 symmetry with a small initial translation.
    // isym1={0,1} and isym2={0,1}: both endpoints are symmetric replicas.
    protein.get_body(0).symmetry().add(symmetry::type::c2);
    protein.get_body(1).symmetry().add(symmetry::type::c2);
    auto* sym_b0 = static_cast<symmetry::CyclicSymmetry*>(protein.get_body(0).symmetry().get(0));
    auto* sym_b1 = static_cast<symmetry::CyclicSymmetry*>(protein.get_body(1).symmetry().get(0));
    sym_b0->_initial_relation.translation = {0.5, 0, 0};
    sym_b1->_initial_relation.translation = {0.5, 0, 0};

    SECTION("relaxed") {
        constraints::DistanceConstraintAtom c(&protein, 0, 0, 1, 0, {0, 1}, {0, 1});
        CHECK(c.evaluate() == 0);
    }

    SECTION("translating body0 changes result") {
        constraints::DistanceConstraintAtom c(&protein, 0, 0, 1, 0, {0, 1}, {0, 1});
        protein.get_body(0).translate(Vector3<double>(1, 0, 0));
        CHECK(c.evaluate() != 0);
        protein.get_body(0).translate(Vector3<double>(-1, 0, 0));
        CHECK(c.evaluate() == 0);
    }

    SECTION("translating body1 changes result") {
        constraints::DistanceConstraintAtom c(&protein, 0, 0, 1, 0, {0, 1}, {0, 1});
        protein.get_body(1).translate(Vector3<double>(1, 0, 0));
        CHECK(c.evaluate() != 0);
        protein.get_body(1).translate(Vector3<double>(-1, 0, 0));
        CHECK(c.evaluate() == 0);
    }

    SECTION("translating unrelated body does not affect result") {
        constraints::DistanceConstraintAtom c(&protein, 0, 0, 1, 0, {0, 1}, {0, 1});
        protein.get_body(2).translate(Vector3<double>(5, 0, 0));
        CHECK(c.evaluate() == 0);
    }

    SECTION("changing tracked symmetry on body0 changes result") {
        constraints::DistanceConstraintAtom c(&protein, 0, 0, 1, 0, {0, 1}, {0, 1});
        sym_b0->_initial_relation.translation = {1.5, 0, 0};
        CHECK(c.evaluate() != 0);
        sym_b0->_initial_relation.translation = {0.5, 0, 0};
        CHECK(c.evaluate() == 0);
    }

    SECTION("changing tracked symmetry on body1 changes result") {
        constraints::DistanceConstraintAtom c(&protein, 0, 0, 1, 0, {0, 1}, {0, 1});
        sym_b1->_initial_relation.translation = {1.5, 0, 0};
        CHECK(c.evaluate() != 0);
        sym_b1->_initial_relation.translation = {0.5, 0, 0};
        CHECK(c.evaluate() == 0);
    }

    SECTION("changing untracked symmetry on body0 does not affect result") {
        protein.get_body(0).symmetry().add(symmetry::type::c2); // sym index 1, not tracked
        auto* sym1_b0 = static_cast<symmetry::CyclicSymmetry*>(protein.get_body(0).symmetry().get(1));
        sym1_b0->_initial_relation.translation = {0.5, 0, 0};
        constraints::DistanceConstraintAtom c(&protein, 0, 0, 1, 0, {0, 1}, {0, 1});
        sym1_b0->_initial_relation.translation = {3, 0, 0}; // modify sym1 - constraint tracks sym0
        CHECK(c.evaluate() == 0);
    }
}

// A constraint may target one of a body's symmetry copies rather than the body itself. The copies are anchored on the body's centre of mass, so the
// constraint has to place them from that centre - not from the representative atom it happens to hold - and it has to keep doing so once the body has
// been rotated away from the frame its symmetry was defined in.
TEST_CASE("DistanceConstraint: a symmetry replica is placed on the body's centre of mass") {
    settings::molecule::implicit_hydrogens = false;

    // exposes the protected distance evaluation so the constraint can be re-checked after the body moves
    struct Probe : constraints::DistanceConstraintAtom {
        using constraints::DistanceConstraintAtom::DistanceConstraintAtom;
        using constraints::DistanceConstraintAtom::evaluate_distance;
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
