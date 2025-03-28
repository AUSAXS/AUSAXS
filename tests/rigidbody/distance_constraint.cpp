#include "settings/MoleculeSettings.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <rigidbody/constraints/DistanceConstraint.h>
#include <data/Molecule.h>
#include <data/Body.h>
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

TEST_CASE_METHOD(fixture, "DistanceConstraint::DistanceConstraint") {
    settings::molecule::implicit_hydrogens = false;
    settings::molecule::center = false;
    Molecule protein(ap);

    SECTION("Protein*, uint, uint, uint, uint") {
        constraints::DistanceConstraint c(&protein, 2, 0, 1, 1);
        CHECK(c.ibody1 == 2);
        CHECK(c.ibody2 == 0);
        CHECK(c.iatom1 == 1);
        CHECK(c.iatom2 == 1);

        CHECK(c.get_body1() == b3);
        CHECK(c.get_body2() == b1);
        CHECK(c.get_atom1() == a6);
        CHECK(c.get_atom2() == a2);
    }

    SECTION("Protein*, Atom&, Atom&") {
        constraints::DistanceConstraint c(&protein, a1, a4);
        CHECK(c.ibody1 == 0);
        CHECK(c.ibody2 == 1);
        CHECK(c.iatom1 == 0);
        CHECK(c.iatom2 == 1);

        CHECK(c.get_body1() == b1);
        CHECK(c.get_body2() == b2);
        CHECK(c.get_atom1() == a1);
        CHECK(c.get_atom2() == a4);
    }

    SECTION("constructor throws with non-carbon atoms") {
        CHECK_THROWS(constraints::DistanceConstraint(&protein, a1, a8));
    }

    SECTION("constructor throws within same body") {
        CHECK_THROWS(constraints::DistanceConstraint(&protein, a1, a2));
    }
}

TEST_CASE_METHOD(fixture, "DistanceConstraint::evaluate") {
    settings::molecule::implicit_hydrogens = false;
    settings::molecule::center = false;
    Molecule protein = Molecule(ap);

    SECTION("relaxed") {
        constraints::DistanceConstraint c13(&protein, a1, a3);
        constraints::DistanceConstraint c24(&protein, a2, a4);
        constraints::DistanceConstraint c35(&protein, a3, a5);
        constraints::DistanceConstraint c47(&protein, a4, a7);
        CHECK(c13.evaluate() == 0);
        CHECK(c24.evaluate() == 0);
        CHECK(c35.evaluate() == 0);
        CHECK(c47.evaluate() == 0);
    }

    SECTION("stretched") {
        constraints::DistanceConstraint c13(&protein, a1, a3);
        constraints::DistanceConstraint c24(&protein, a2, a4);
        constraints::DistanceConstraint c35(&protein, a3, a5);
        constraints::DistanceConstraint c57(&protein, a5, a7);

        protein.get_body(0).translate(Vector3<double>(1, 0, 0));
        CHECK(c13.evaluate() != 0);
        CHECK(c24.evaluate() != 0);
        CHECK(c35.evaluate() == 0);
        CHECK(c57.evaluate() == 0);

        protein.get_body(0).translate(Vector3<double>(-1, 0, 0));
        CHECK(c13.evaluate() == 0);
        CHECK(c24.evaluate() == 0);
        CHECK(c35.evaluate() == 0);
        CHECK(c57.evaluate() == 0);

        protein.get_body(1).translate(Vector3<double>(1, 0, 0));
        CHECK(c13.evaluate() != 0);
        CHECK(c24.evaluate() != 0);
        CHECK(c35.evaluate() != 0);
        CHECK(c57.evaluate() == 0);

        protein.get_body(3).translate(Vector3<double>(1, 0, 0));
        CHECK(c57.evaluate() != 0);
    }
}

TEST_CASE_METHOD(fixture, "DistanceConstraint::operator==") {
    settings::molecule::implicit_hydrogens = false;
    settings::molecule::center = false;
    Molecule protein = Molecule(ap);

    constraints::DistanceConstraint c1(&protein, a1, a3);
    constraints::DistanceConstraint c2(&protein, a1, a3);
    constraints::DistanceConstraint c3(&protein, a1, a4);

    CHECK(c1 == c2);
    CHECK_FALSE(c1 == c3);
}