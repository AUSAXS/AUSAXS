#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <rigidbody/constraints/OverlapConstraint.h>
#include <rigidbody/constraints/DistanceConstraint.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <data/Protein.h>
#include <data/Body.h>
#include <data/Atom.h>
#include <settings/All.h>

using namespace rigidbody;

struct fixture {
    fixture() {
        settings::protein::use_effective_charge = false;
        settings::axes::distance_bin_width = 0.1;
        settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::None;
    }

    Atom a1 = Atom(Vector3<double>(-1, -1, -1), 1, "C", "C", 1);
    Atom a2 = Atom(Vector3<double>(-1,  1, -1), 1, "C", "C", 1);
    Atom a3 = Atom(Vector3<double>(-1, -1,  1), 1, "C", "C", 1);
    Atom a4 = Atom(Vector3<double>(-1,  1,  1), 1, "C", "C", 1);
    Atom a5 = Atom(Vector3<double>( 1, -1, -1), 1, "C", "C", 1);
    Atom a6 = Atom(Vector3<double>( 1,  1, -1), 1, "C", "C", 1);
    Atom a7 = Atom(Vector3<double>( 1, -1,  1), 1, "C", "C", 1);
    Atom a8 = Atom(Vector3<double>( 1,  1,  1), 1, "He", "He", 1);

    Body b1 = Body(std::vector<Atom>{a1, a2});
    Body b2 = Body(std::vector<Atom>{a3, a4});
    Body b3 = Body(std::vector<Atom>{a5, a6});
    Body b4 = Body(std::vector<Atom>{a7, a8});
    std::vector<Body> ap = {b1, b2, b3, b4};
};

TEST_CASE_METHOD(fixture, "ConstraintManager::ConstraintManager") {
    Protein protein(ap);
    SECTION("Protein*") {
        ConstraintManager cm(&protein);
        CHECK(cm.protein == &protein);
    }
}

TEST_CASE_METHOD(fixture, "ConstraintManager::add_constraint") {
    Protein protein(ap);
    SECTION("OverlapConstraint&&") {
        ConstraintManager cm(&protein);
        OverlapConstraint oc(&protein);
        auto oc_copy = oc;
        cm.add_constraint(std::move(oc));
        CHECK(cm.overlap_constraint == oc_copy);
    }

    SECTION("OverlapConstraint&") {
        ConstraintManager cm(&protein);
        OverlapConstraint oc(&protein);
        cm.add_constraint(oc);
        CHECK(cm.overlap_constraint == oc);
    }

    SECTION("DistanceConstraint&&") {
        ConstraintManager cm(&protein);
        DistanceConstraint dc(&protein, a1, a3);
        auto dc_copy = dc;
        cm.add_constraint(std::move(dc));
        CHECK(cm.distance_constraints.size() == 1);
        CHECK(cm.distance_constraints[0] == dc_copy);
    }

    SECTION("DistanceConstraint&") {
        ConstraintManager cm(&protein);
        DistanceConstraint dc(&protein, a1, a3);
        cm.add_constraint(dc);
        CHECK(cm.distance_constraints.size() == 1);
        CHECK(cm.distance_constraints[0] == dc);
    }

    SECTION("Multiple") {
        ConstraintManager cm(&protein);
        OverlapConstraint oc(&protein);
        DistanceConstraint dc1(&protein, a1, a3);
        DistanceConstraint dc2(&protein, a1, a4);
        cm.add_constraint(oc);
        cm.add_constraint(dc1);
        cm.add_constraint(std::move(dc2));
        CHECK(cm.overlap_constraint == oc);
        REQUIRE(cm.distance_constraints.size() == 2);
        CHECK(cm.distance_constraints[0] == dc1);
        CHECK(cm.distance_constraints[1] == dc2);
    }
}

TEST_CASE_METHOD(fixture, "ConstraintManager::evaluate") {
    Protein protein(ap);
    SECTION("returns chi2 contribution of all constraints") {
        ConstraintManager cm(&protein);
        CHECK(cm.evaluate() == 0);

        OverlapConstraint oc(&protein);
        DistanceConstraint dc1(&protein, a1, a3);
        DistanceConstraint dc2(&protein, a1, a4);

        cm.add_constraint(oc);
        protein.get_body(0).translate(Vector3<double>(2, 2, 1.5));
        CHECK(oc.evaluate() != 0);
        CHECK_THAT(cm.evaluate(), Catch::Matchers::WithinAbs(oc.evaluate(), 1e-3));
        protein.get_body(0).translate(Vector3<double>(-2, -2, -1.5));

        cm.add_constraint(dc1);
        protein.get_body(0).translate(Vector3<double>(2, 2, 1.5));
        CHECK(oc.evaluate() != 0);
        CHECK(dc1.evaluate() != 0);
        CHECK_THAT(cm.evaluate(), Catch::Matchers::WithinAbs(oc.evaluate() + dc1.evaluate(), 1e-3));
        protein.get_body(0).translate(Vector3<double>(-2, -2, -1.5));

        cm.add_constraint(std::move(dc2));
        protein.get_body(0).translate(Vector3<double>(2, 2, 1.5));
        CHECK(oc.evaluate() != 0);
        CHECK(dc1.evaluate() != 0);
        CHECK(dc2.evaluate() != 0);
        CHECK_THAT(cm.evaluate(), Catch::Matchers::WithinAbs(oc.evaluate() + dc1.evaluate() + dc2.evaluate(), 1e-3));
    }
}