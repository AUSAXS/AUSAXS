#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <rigidbody/constraints/DistanceConstraintAtom.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <data/symmetry/PredefinedSymmetries.h>
#include <data/symmetry/CyclicSymmetry.h>
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
    settings::molecule::center = false;
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
    settings::molecule::center = false;
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
}

TEST_CASE_METHOD(fixture, "DistanceConstraintAtom::evaluate with symmetry") {
    settings::molecule::implicit_hydrogens = false;
    settings::molecule::center = false;
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
}