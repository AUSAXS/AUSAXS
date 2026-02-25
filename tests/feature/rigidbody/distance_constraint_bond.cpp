// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <rigidbody/constraints/DistanceConstraintBond.h>
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
    AtomFF a8 = AtomFF({ 1,  1,  1}, form_factor::form_factor_t::NH);

    Body b1 = Body(std::vector<AtomFF>{a1, a2});
    Body b2 = Body(std::vector<AtomFF>{a3, a4});
    Body b3 = Body(std::vector<AtomFF>{a5, a8}); // b3 has only one C atom
    std::vector<Body> ap = {b1, b2, b3};
};

TEST_CASE_METHOD(fixture, "DistanceConstraintBond::constructor") {
    settings::molecule::implicit_hydrogens = false;
    settings::molecule::center = false;
    Molecule protein = Molecule(ap);

    SECTION("selects the closest C-C pair between the two bodies") {
        constraints::DistanceConstraintBond c(&protein, 0, 1);
        // a1(-1,-1,-1) to a3(-1,-1,1): dist=2 is the closest pair (found first)
        CHECK(c.ibody1 == 0);
        CHECK(c.ibody2 == 1);
        CHECK(c.get_atom1() == a1);
        CHECK(c.get_atom2() == a3);
        CHECK_THAT(c.d_target, Catch::Matchers::WithinAbs(2.0, 1e-9));
    }

    SECTION("throws if atoms are too far apart") {
        Body far1 = Body(std::vector<AtomFF>{AtomFF({0,  0, 0}, form_factor::form_factor_t::C)});
        Body far2 = Body(std::vector<AtomFF>{AtomFF({0, 10, 0}, form_factor::form_factor_t::C)});
        std::vector<Body> far_ap = {far1, far2};
        Molecule far_protein(far_ap);
        CHECK_THROWS(constraints::DistanceConstraintBond(&far_protein, 0, 1));
    }

    SECTION("stores symmetry indices") {
        protein.get_body(0).symmetry().add(symmetry::type::c2);
        constraints::DistanceConstraintBond c(&protein, 0, 1, {0, 1}, {-1, -1});
        CHECK(c.isym1 == std::make_pair(0, 1));
        CHECK(c.isym2 == std::make_pair(-1, -1));
    }
}

TEST_CASE_METHOD(fixture, "DistanceConstraintBond::evaluate") {
    settings::molecule::implicit_hydrogens = false;
    settings::molecule::center = false;
    Molecule protein = Molecule(ap);

    SECTION("relaxed") {
        constraints::DistanceConstraintBond c(&protein, 0, 1);
        CHECK(c.evaluate() == 0);
    }

    SECTION("stretched - penalty above d_target") {
        constraints::DistanceConstraintBond c(&protein, 0, 1);

        // Move body0 away from body1 along z: distance increases beyond d_target
        protein.get_body(0).translate(Vector3<double>(0, 0, -2));
        CHECK(c.evaluate() != 0);
    }

    SECTION("compressed - penalty below d_target") {
        constraints::DistanceConstraintBond c(&protein, 0, 1);

        // Move body0 towards body1 along z: distance decreases below d_target
        protein.get_body(0).translate(Vector3<double>(0, 0, 1));
        CHECK(c.evaluate() != 0);

        protein.get_body(0).translate(Vector3<double>(0, 0, -1));
        CHECK(c.evaluate() == 0);
    }

    SECTION("unrelated body translation does not affect constraint") {
        constraints::DistanceConstraintBond c(&protein, 0, 1);
        protein.get_body(2).translate(Vector3<double>(5, 0, 0));
        CHECK(c.evaluate() == 0);
    }

    SECTION("symmetries on unrelated body do not affect constraint") {
        constraints::DistanceConstraintBond c(&protein, 0, 1);
        protein.get_body(2).symmetry().add(symmetry::type::c2);
        auto* sym2 = static_cast<symmetry::CyclicSymmetry*>(protein.get_body(2).symmetry().get(0));
        sym2->_initial_relation.translation = {0.5, 0, 0};
        CHECK(c.evaluate() == 0);
    }
}

TEST_CASE_METHOD(fixture, "DistanceConstraintBond::evaluate with symmetry") {
    settings::molecule::implicit_hydrogens = false;
    settings::molecule::center = false;
    Molecule protein = Molecule(ap);

    // Use a small initial translation so the symmetric replica of body0's atoms
    // stays within 4Å of body1's atoms (Bond constructor throws otherwise).
    protein.get_body(0).symmetry().add(symmetry::type::c2);
    auto* sym = static_cast<symmetry::CyclicSymmetry*>(protein.get_body(0).symmetry().get(0));
    sym->_initial_relation.translation = {0.5, 0, 0};
    // Symmetric copy of a1(-1,-1,-1) lands at (-2,-1,-1), which is ~2.24Å from a3(-1,-1,1)

    SECTION("relaxed") {
        constraints::DistanceConstraintBond c(&protein, 0, 1, {0, 1}, {-1, -1});
        CHECK(c.evaluate() == 0);
    }

    SECTION("stretched - penalty") {
        constraints::DistanceConstraintBond c(&protein, 0, 1, {0, 1}, {-1, -1});

        // Move body0 away from body1 along z: symmetric copy also moves away
        protein.get_body(0).translate(Vector3<double>(0, 0, -1));
        CHECK(c.evaluate() != 0);
    }

    SECTION("compressed") {
        constraints::DistanceConstraintBond c(&protein, 0, 1, {0, 1}, {-1, -1});

        // Move body0 toward body1 along z: symmetric copy gets closer than d_target
        protein.get_body(0).translate(Vector3<double>(0, 0, 1));
        CHECK(c.evaluate() != 0);

        protein.get_body(0).translate(Vector3<double>(0, 0, -1));
        CHECK(c.evaluate() == 0);
    }

    SECTION("changing tracked symmetry translation changes result") {
        constraints::DistanceConstraintBond c(&protein, 0, 1, {0, 1}, {-1, -1});
        sym->_initial_relation.translation = {2.0, 0, 0}; // was 0.5
        CHECK(c.evaluate() != 0);
        sym->_initial_relation.translation = {0.5, 0, 0};
        CHECK(c.evaluate() == 0);
    }

    SECTION("changing untracked symmetry on body0 does not affect result") {
        protein.get_body(0).symmetry().add(symmetry::type::c2); // sym index 1, not tracked
        auto* sym1 = static_cast<symmetry::CyclicSymmetry*>(protein.get_body(0).symmetry().get(1));
        sym1->_initial_relation.translation = {0.5, 0, 0};
        constraints::DistanceConstraintBond c(&protein, 0, 1, {0, 1}, {-1, -1});
        sym1->_initial_relation.translation = {2.0, 0, 0}; // modify sym1 - constraint tracks sym0
        CHECK(c.evaluate() == 0);
    }

    SECTION("symmetry and translation on unrelated body does not affect result") {
        constraints::DistanceConstraintBond c(&protein, 0, 1, {0, 1}, {-1, -1});
        protein.get_body(2).symmetry().add(symmetry::type::c2);
        auto* sym2 = static_cast<symmetry::CyclicSymmetry*>(protein.get_body(2).symmetry().get(0));
        sym2->_initial_relation.translation = {0.5, 0, 0};
        protein.get_body(2).translate(Vector3<double>(5, 0, 0));
        CHECK(c.evaluate() == 0);
    }
}

TEST_CASE_METHOD(fixture, "DistanceConstraintBond::evaluate symmetry-symmetry") {
    settings::molecule::implicit_hydrogens = false;
    settings::molecule::center = false;
    Molecule protein = Molecule(ap);

    // Both body0 and body1 get a c2 symmetry with a small initial translation.
    // With t_i={0.5,0,0}: transformed a1(-1,-1,-1)->(-2,-1,-1), a3(-1,-1,1)->(-2,-1,1), d=2 ≤ 4.
    protein.get_body(0).symmetry().add(symmetry::type::c2);
    protein.get_body(1).symmetry().add(symmetry::type::c2);
    auto* sym_b0 = static_cast<symmetry::CyclicSymmetry*>(protein.get_body(0).symmetry().get(0));
    auto* sym_b1 = static_cast<symmetry::CyclicSymmetry*>(protein.get_body(1).symmetry().get(0));
    sym_b0->_initial_relation.translation = {0.5, 0, 0};
    sym_b1->_initial_relation.translation = {0.5, 0, 0};

    SECTION("relaxed") {
        constraints::DistanceConstraintBond c(&protein, 0, 1, {0, 1}, {0, 1});
        CHECK(c.evaluate() == 0);
    }

    SECTION("translating body0 changes result") {
        constraints::DistanceConstraintBond c(&protein, 0, 1, {0, 1}, {0, 1});
        protein.get_body(0).translate(Vector3<double>(0, 0, 1));
        CHECK(c.evaluate() != 0);
        protein.get_body(0).translate(Vector3<double>(0, 0, -1));
        CHECK(c.evaluate() == 0);
    }

    SECTION("translating body1 changes result") {
        constraints::DistanceConstraintBond c(&protein, 0, 1, {0, 1}, {0, 1});
        protein.get_body(1).translate(Vector3<double>(0, 0, 1));
        CHECK(c.evaluate() != 0);
        protein.get_body(1).translate(Vector3<double>(0, 0, -1));
        CHECK(c.evaluate() == 0);
    }

    SECTION("translating unrelated body does not affect result") {
        constraints::DistanceConstraintBond c(&protein, 0, 1, {0, 1}, {0, 1});
        protein.get_body(2).translate(Vector3<double>(5, 0, 0));
        CHECK(c.evaluate() == 0);
    }

    SECTION("changing tracked symmetry on body0 changes result") {
        constraints::DistanceConstraintBond c(&protein, 0, 1, {0, 1}, {0, 1});
        sym_b0->_initial_relation.translation = {1.5, 0, 0};
        CHECK(c.evaluate() != 0);
        sym_b0->_initial_relation.translation = {0.5, 0, 0};
        CHECK(c.evaluate() == 0);
    }

    SECTION("changing tracked symmetry on body1 changes result") {
        constraints::DistanceConstraintBond c(&protein, 0, 1, {0, 1}, {0, 1});
        sym_b1->_initial_relation.translation = {1.5, 0, 0};
        CHECK(c.evaluate() != 0);
        sym_b1->_initial_relation.translation = {0.5, 0, 0};
        CHECK(c.evaluate() == 0);
    }

    SECTION("changing untracked symmetry on body0 does not affect result") {
        protein.get_body(0).symmetry().add(symmetry::type::c2); // sym index 1, not tracked
        auto* sym1_b0 = static_cast<symmetry::CyclicSymmetry*>(protein.get_body(0).symmetry().get(1));
        sym1_b0->_initial_relation.translation = {0.5, 0, 0};
        constraints::DistanceConstraintBond c(&protein, 0, 1, {0, 1}, {0, 1});
        sym1_b0->_initial_relation.translation = {3, 0, 0}; // modify sym1 - constraint tracks sym0
        CHECK(c.evaluate() == 0);
    }
}
