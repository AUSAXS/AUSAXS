#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <rigidbody/constraints/DistanceConstraintCM.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <data/symmetry/PredefinedSymmetries.h>
#include <data/symmetry/CyclicSymmetry.h>
#include <form_factor/FormFactor.h>
#include <constants/Constants.h>
#include <settings/All.h>

using namespace ausaxs;
using namespace ausaxs::data;
using namespace ausaxs::rigidbody;

struct fixture {
    // Single-C bodies: only one valid atom per body, selection is unambiguous.
    AtomFF c1  = AtomFF({0, 0, 0}, form_factor::form_factor_t::C);
    AtomFF nh1 = AtomFF({8, 0, 0}, form_factor::form_factor_t::NH);
    AtomFF c2  = AtomFF({0, 0, 4}, form_factor::form_factor_t::C);
    AtomFF nh2 = AtomFF({8, 0, 4}, form_factor::form_factor_t::NH);

    Body b1 = Body(std::vector<AtomFF>{c1, nh1});
    Body b2 = Body(std::vector<AtomFF>{c2, nh2});
    std::vector<Body> ap = {b1, b2};
};

TEST_CASE_METHOD(fixture, "DistanceConstraintCM::constructor") {
    settings::molecule::implicit_hydrogens = false;
    settings::molecule::center = false;
    Molecule protein = Molecule(ap);

    SECTION("selects only C atoms, not NH") {
        constraints::DistanceConstraintCM c(&protein, 0, 1);
        CHECK(c.ibody1 == 0);
        CHECK(c.ibody2 == 1);
        CHECK(form_factor::to_atom_type(c.get_atom1().form_factor_type()) == constants::atom_t::C);
        CHECK(form_factor::to_atom_type(c.get_atom2().form_factor_type()) == constants::atom_t::C);
        // With a single C per body, the unique C must be selected
        CHECK(c.get_atom1() == c1);
        CHECK(c.get_atom2() == c2);
    }

    SECTION("d_target is set to current distance at construction") {
        constraints::DistanceConstraintCM c(&protein, 0, 1);
        CHECK_THAT(c.d_target, Catch::Matchers::WithinAbs(4.0, 1e-9)); // distance c1(0,0,0) to c2(0,0,4)
    }

    SECTION("selects C atom closest to CM in multi-C body") {
        // Three C atoms: two at extremes and one near CM.
        // Regardless of mass, the middle atom (near the average position) is closest to CM.
        AtomFF c_left  = AtomFF({-10, 0, 0}, form_factor::form_factor_t::C);
        AtomFF c_mid   = AtomFF({  0, 0, 0}, form_factor::form_factor_t::C);
        AtomFF c_right = AtomFF({  1, 0, 0}, form_factor::form_factor_t::C);
        Body b1_mc = Body(std::vector<AtomFF>{c_left, c_mid, c_right});
        Body b2_s  = Body(std::vector<AtomFF>{AtomFF({0, 5, 0}, form_factor::form_factor_t::C)});
        Molecule protein2(std::vector<Body>{b1_mc, b2_s});

        constraints::DistanceConstraintCM c2(&protein2, 0, 1);
        // CM of b1_mc is near (-10+0+1)/3 = -3; c_mid at 0 is at distance 3, c_right at distance 4, c_left at distance 7
        CHECK(c2.get_atom1() == c_mid);
    }

    SECTION("stores symmetry indices") {
        protein.get_body(0).symmetry().add(symmetry::type::c2);
        constraints::DistanceConstraintCM c(&protein, 0, 1, {0, 1}, {-1, -1});
        CHECK(c.isym1 == std::make_pair(0, 1));
        CHECK(c.isym2 == std::make_pair(-1, -1));
    }
}

TEST_CASE_METHOD(fixture, "DistanceConstraintCM::evaluate") {
    settings::molecule::implicit_hydrogens = false;
    settings::molecule::center = false;
    Molecule protein = Molecule(ap);

    SECTION("relaxed") {
        constraints::DistanceConstraintCM c(&protein, 0, 1);
        CHECK(c.evaluate() == 0);
    }

    SECTION("stretched") {
        constraints::DistanceConstraintCM c(&protein, 0, 1);

        protein.get_body(0).translate(Vector3<double>(0, 0, -2));
        CHECK(c.evaluate() != 0);

        protein.get_body(0).translate(Vector3<double>(0, 0, 2));
        CHECK(c.evaluate() == 0);
    }

    SECTION("compressed") {
        constraints::DistanceConstraintCM c(&protein, 0, 1);

        protein.get_body(0).translate(Vector3<double>(0, 0, 2));
        CHECK(c.evaluate() != 0);

        protein.get_body(0).translate(Vector3<double>(0, 0, -2));
        CHECK(c.evaluate() == 0);
    }

    SECTION("unrelated body translation does not affect constraint") {
        Body b3 = Body(std::vector<AtomFF>{AtomFF({10, 0, 0}, form_factor::form_factor_t::C)});
        std::vector<Body> ap3 = {b1, b2, b3};
        Molecule protein2(ap3);
        constraints::DistanceConstraintCM c(&protein2, 0, 1);

        protein2.get_body(2).translate(Vector3<double>(5, 0, 0));
        CHECK(c.evaluate() == 0);
    }

    SECTION("symmetries on unrelated body do not affect real-real constraint") {
        Body b3 = Body(std::vector<AtomFF>{AtomFF({10, 0, 0}, form_factor::form_factor_t::C)});
        std::vector<Body> ap3 = {b1, b2, b3};
        Molecule protein2(ap3);
        constraints::DistanceConstraintCM c(&protein2, 0, 1);

        protein2.get_body(2).symmetry().add(symmetry::type::c2);
        auto* sym2 = static_cast<symmetry::CyclicSymmetry*>(protein2.get_body(2).symmetry().get(0));
        sym2->_initial_relation.translation = {3, 0, 0};
        CHECK(c.evaluate() == 0);

        protein2.get_body(2).translate(Vector3<double>(5, 0, 0));
        CHECK(c.evaluate() == 0);
    }
}

TEST_CASE_METHOD(fixture, "DistanceConstraintCM::evaluate with symmetry") {
    settings::molecule::implicit_hydrogens = false;
    settings::molecule::center = false;
    Molecule protein = Molecule(ap);

    protein.get_body(0).symmetry().add(symmetry::type::c2);
    auto* sym = static_cast<symmetry::CyclicSymmetry*>(protein.get_body(0).symmetry().get(0));
    sym->_initial_relation.translation = {5, 0, 0};

    SECTION("relaxed") {
        constraints::DistanceConstraintCM c(&protein, 0, 1, {0, 1}, {-1, -1});
        CHECK(c.evaluate() == 0);
    }

    SECTION("stretched") {
        constraints::DistanceConstraintCM c(&protein, 0, 1, {0, 1}, {-1, -1});

        protein.get_body(0).translate(Vector3<double>(1, 0, 0));
        CHECK(c.evaluate() != 0);

        protein.get_body(0).translate(Vector3<double>(-1, 0, 0));
        CHECK(c.evaluate() == 0);
    }

    SECTION("changing tracked symmetry translation changes result") {
        constraints::DistanceConstraintCM c(&protein, 0, 1, {0, 1}, {-1, -1});
        sym->_initial_relation.translation = {8, 0, 0}; // was 5
        CHECK(c.evaluate() != 0);
        sym->_initial_relation.translation = {5, 0, 0};
        CHECK(c.evaluate() == 0);
    }

    SECTION("changing untracked symmetry on constrained body does not affect result") {
        protein.get_body(0).symmetry().add(symmetry::type::c2); // sym index 1, not tracked
        auto* sym1 = static_cast<symmetry::CyclicSymmetry*>(protein.get_body(0).symmetry().get(1));
        sym1->_initial_relation.translation = {0.5, 0, 0};
        constraints::DistanceConstraintCM c(&protein, 0, 1, {0, 1}, {-1, -1});
        sym1->_initial_relation.translation = {9, 0, 0}; // modify sym1 - constraint tracks sym0
        CHECK(c.evaluate() == 0);
    }

    SECTION("symmetry and translation on unrelated body does not affect result") {
        Body b3 = Body(std::vector<AtomFF>{AtomFF({10, 0, 0}, form_factor::form_factor_t::C)});
        Molecule protein2(std::vector<Body>{b1, b2, b3});
        protein2.get_body(0).symmetry().add(symmetry::type::c2);
        auto* sym2 = static_cast<symmetry::CyclicSymmetry*>(protein2.get_body(0).symmetry().get(0));
        sym2->_initial_relation.translation = {5, 0, 0};
        constraints::DistanceConstraintCM c(&protein2, 0, 1, {0, 1}, {-1, -1});

        protein2.get_body(2).symmetry().add(symmetry::type::c2);
        auto* sym3 = static_cast<symmetry::CyclicSymmetry*>(protein2.get_body(2).symmetry().get(0));
        sym3->_initial_relation.translation = {3, 0, 0};
        protein2.get_body(2).translate(Vector3<double>(5, 0, 0));
        CHECK(c.evaluate() == 0);
    }
}

TEST_CASE_METHOD(fixture, "DistanceConstraintCM::evaluate symmetry-symmetry") {
    settings::molecule::implicit_hydrogens = false;
    settings::molecule::center = false;
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
        constraints::DistanceConstraintCM c(&protein, 0, 1, {0, 1}, {0, 1});
        CHECK(c.evaluate() == 0);
    }

    SECTION("translating body0 changes result") {
        constraints::DistanceConstraintCM c(&protein, 0, 1, {0, 1}, {0, 1});
        protein.get_body(0).translate(Vector3<double>(1, 0, 0));
        CHECK(c.evaluate() != 0);
        protein.get_body(0).translate(Vector3<double>(-1, 0, 0));
        CHECK(c.evaluate() == 0);
    }

    SECTION("translating body1 changes result") {
        constraints::DistanceConstraintCM c(&protein, 0, 1, {0, 1}, {0, 1});
        protein.get_body(1).translate(Vector3<double>(1, 0, 0));
        CHECK(c.evaluate() != 0);
        protein.get_body(1).translate(Vector3<double>(-1, 0, 0));
        CHECK(c.evaluate() == 0);
    }

    SECTION("translating unrelated body does not affect result") {
        Body b3 = Body(std::vector<AtomFF>{AtomFF({10, 0, 0}, form_factor::form_factor_t::C)});
        Molecule protein2(std::vector<Body>{b1, b2, b3});
        protein2.get_body(0).symmetry().add(symmetry::type::c2);
        protein2.get_body(1).symmetry().add(symmetry::type::c2);
        auto* sym2_b0 = static_cast<symmetry::CyclicSymmetry*>(protein2.get_body(0).symmetry().get(0));
        auto* sym2_b1 = static_cast<symmetry::CyclicSymmetry*>(protein2.get_body(1).symmetry().get(0));
        sym2_b0->_initial_relation.translation = {0.5, 0, 0};
        sym2_b1->_initial_relation.translation = {0.5, 0, 0};
        constraints::DistanceConstraintCM c(&protein2, 0, 1, {0, 1}, {0, 1});

        protein2.get_body(2).translate(Vector3<double>(5, 0, 0));
        CHECK(c.evaluate() == 0);
    }

    SECTION("changing tracked symmetry on body0 changes result") {
        constraints::DistanceConstraintCM c(&protein, 0, 1, {0, 1}, {0, 1});
        sym_b0->_initial_relation.translation = {1.5, 0, 0};
        CHECK(c.evaluate() != 0);
        sym_b0->_initial_relation.translation = {0.5, 0, 0};
        CHECK(c.evaluate() == 0);
    }

    SECTION("changing tracked symmetry on body1 changes result") {
        constraints::DistanceConstraintCM c(&protein, 0, 1, {0, 1}, {0, 1});
        sym_b1->_initial_relation.translation = {1.5, 0, 0};
        CHECK(c.evaluate() != 0);
        sym_b1->_initial_relation.translation = {0.5, 0, 0};
        CHECK(c.evaluate() == 0);
    }

    SECTION("changing untracked symmetry on body0 does not affect result") {
        protein.get_body(0).symmetry().add(symmetry::type::c2); // sym index 1, not tracked
        auto* sym1_b0 = static_cast<symmetry::CyclicSymmetry*>(protein.get_body(0).symmetry().get(1));
        sym1_b0->_initial_relation.translation = {0.5, 0, 0};
        constraints::DistanceConstraintCM c(&protein, 0, 1, {0, 1}, {0, 1});
        sym1_b0->_initial_relation.translation = {3, 0, 0}; // modify sym1 - constraint tracks sym0
        CHECK(c.evaluate() == 0);
    }
}
