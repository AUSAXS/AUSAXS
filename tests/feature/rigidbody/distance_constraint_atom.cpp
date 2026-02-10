#include "settings/MoleculeSettings.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <rigidbody/constraints/DistanceConstraintAtom.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <data/symmetry/PredefinedSymmetries.h>
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

    SECTION("constructor with symmetry stores indices") {
        // Use the 2-body constructor which does not evaluate distance, to test index storage
        constraints::DistanceConstraintAtom c(&protein, 0, 1, {0, -1}, {1, -1});
        CHECK(c.ibody1 == 0);
        CHECK(c.ibody2 == 1);
        CHECK(c.isym1 == std::make_pair(0, -1));
        CHECK(c.isym2 == std::make_pair(1, -1));
    }

    SECTION("constructor throws with non-carbon atoms") {
        CHECK_THROWS(constraints::DistanceConstraintAtom(&protein, 0, 0, 3, 1));
    }

    SECTION("constructor throws within same body") {
        CHECK_THROWS(constraints::DistanceConstraintAtom(&protein, 0, 0, 0, 1));
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
}

TEST_CASE_METHOD(fixture, "DistanceConstraintAtom::evaluate with symmetry", "[broken]") {
    settings::molecule::implicit_hydrogens = false;
    settings::molecule::center = false;
    Molecule protein = Molecule(ap);

    // Add dimer symmetry to body 0
    protein.get_body(0).symmetry().add(symmetry::type::p2);

    SECTION("symmetric constraint") {
        // Same-body symmetry constraints not supported by current DistanceConstraintAtom constructor
    }

    SECTION("symmetric constraint with translation") {
        // Same-body symmetry constraints not supported by current DistanceConstraintAtom constructor
    }
}