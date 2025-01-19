#include "settings/MoleculeSettings.h"
#include <catch2/catch_test_macros.hpp>

#include <rigidbody/constraints/DistanceConstraint.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/transform/RigidTransform.h>
#include <rigidbody/transform/SingleTransform.h>
#include <rigidbody/transform/TransformGroup.h>

#include <data/Body.h>
#include <rigidbody/RigidBody.h>
#include <math/MatrixUtils.h>
#include <settings/All.h>
#include <unordered_set>

using namespace ausaxs;
using namespace ausaxs::data;
using namespace ausaxs::rigidbody;

struct fixture {
    fixture() {
        settings::molecule::implicit_hydrogens = false;
        settings::molecule::center = false;

        a1 = AtomFF({-1, -1, -1}, form_factor::form_factor_t::C);
        a2 = AtomFF({-1,  1, -1}, form_factor::form_factor_t::C);
        b1 = Body(std::vector{a1, a2});

        a3 = AtomFF({1, -1, -1}, form_factor::form_factor_t::C);
        a4 = AtomFF({1,  1, -1}, form_factor::form_factor_t::C);
        b2 = Body(std::vector{a3, a4});

        a5 =  AtomFF({-1, -1,  1}, form_factor::form_factor_t::C);
        a6 =  AtomFF({-1,  1,  1}, form_factor::form_factor_t::C);
        b3 = Body(std::vector{a5, a6});

        a7 =  AtomFF({ 1, -1,  1}, form_factor::form_factor_t::C);
        a8 =  AtomFF({ 1,  1,  1}, form_factor::form_factor_t::C);
        b4 = Body(std::vector{a7, a8});

        a9 =  AtomFF({ 0,  0,  0}, form_factor::form_factor_t::C);
        a10 = AtomFF({ 0,  0,  2}, form_factor::form_factor_t::C);
        b5 = Body(std::vector{a9, a10});

        bodies = {b1, b2, b3, b4, b5};
    }

    AtomFF a1, a2, a3, a4, a5, a6, a7, a8, a9, a10;
    Body b1, b2, b3, b4, b5;
    std::vector<Body> bodies;
};

TEST_CASE_METHOD(fixture, "TransformStrategy::apply") {
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::None;
    settings::general::verbose = false;
    settings::grid::scaling = 2;

    SECTION("SingleTransform::apply") {
        RigidBody rigidbody(bodies);
        auto manager = rigidbody.get_constraint_manager();

        manager->add_constraint(rigidbody::constraints::DistanceConstraint(&rigidbody, 0, 1, 0, 0)); // 0 <-- this one
        manager->add_constraint(rigidbody::constraints::DistanceConstraint(&rigidbody, 1, 2, 0, 0)); // 1
        manager->add_constraint(rigidbody::constraints::DistanceConstraint(&rigidbody, 2, 3, 0, 0)); // 2
        manager->add_constraint(rigidbody::constraints::DistanceConstraint(&rigidbody, 3, 4, 0, 0)); // 3

        // translate
        rigidbody::transform::SingleTransform transform(&rigidbody);
        transform.apply({Matrix<double>::identity(3), Vector3<double>(1, 0, 0)}, manager->distance_constraints[0]);
        CHECK(rigidbody.get_body(0).get_atom(0).coordinates() == Vector3<double>(0, -1, -1));
        CHECK(rigidbody.get_body(0).get_atom(1).coordinates() == Vector3<double>(0,  1, -1));
        transform.undo();
        CHECK(rigidbody.get_body(0).get_atom(0).coordinates() == Vector3<double>(-1, -1, -1));
        CHECK(rigidbody.get_body(0).get_atom(1).coordinates() == Vector3<double>(-1,  1, -1));

        // rotate
        transform.apply({matrix::rotation_matrix(Vector3<double>{0, 0, 1}, M_PI/2), Vector3<double>(0, 0, 0)}, manager->distance_constraints[0]);
        CHECK(rigidbody.get_body(0).get_atom(0).coordinates() == Vector3<double>( 1, -3, -1));
        CHECK(rigidbody.get_body(0).get_atom(1).coordinates() == Vector3<double>(-1, -3, -1));
        transform.undo();
        CHECK(rigidbody.get_body(0).get_atom(0).coordinates() == Vector3<double>(-1, -1, -1));
        CHECK(rigidbody.get_body(0).get_atom(1).coordinates() == Vector3<double>(-1,  1, -1));
    }

    SECTION("RigidTransform::apply") {
        // generate an easier constraint to test the apply function on
        AtomFF a11({3,  1, -1}, form_factor::form_factor_t::C);
        AtomFF a12({3, -1,  1}, form_factor::form_factor_t::C);
        Body b6(std::vector<AtomFF>{a11, a12});
        RigidBody rigidbody(std::vector<Body>{b1, b2, b3, b4, b5, b6});
        auto manager = rigidbody.get_constraint_manager();

        manager->add_constraint(rigidbody::constraints::DistanceConstraint(&rigidbody, 0, 1, 0, 0)); // 0 <-- first this one
        manager->add_constraint(rigidbody::constraints::DistanceConstraint(&rigidbody, 1, 5, 1, 0)); // 1 <-- then this one
        manager->add_constraint(rigidbody::constraints::DistanceConstraint(&rigidbody, 5, 2, 1, 0)); // 2
        manager->add_constraint(rigidbody::constraints::DistanceConstraint(&rigidbody, 2, 3, 0, 0)); // 3
        manager->add_constraint(rigidbody::constraints::DistanceConstraint(&rigidbody, 3, 4, 0, 0)); // 4

        // single-body translate
        rigidbody::transform::RigidTransform transform(&rigidbody);
        transform.apply({Matrix<double>::identity(3), Vector3<double>(1, 0, 0)}, manager->distance_constraints[0]);
        CHECK(rigidbody.get_body(0).get_atom(0).coordinates() == Vector3<double>(0, -1, -1));
        CHECK(rigidbody.get_body(0).get_atom(1).coordinates() == Vector3<double>(0,  1, -1));
        transform.undo();
        CHECK(rigidbody.get_body(0).get_atom(0).coordinates() == Vector3<double>(-1, -1, -1));
        CHECK(rigidbody.get_body(0).get_atom(1).coordinates() == Vector3<double>(-1,  1, -1));

        // single-body rotate
        transform.apply({matrix::rotation_matrix(Vector3<double>{0, 0, 1}, M_PI/2), Vector3<double>(0, 0, 0)}, manager->distance_constraints[0]);
        CHECK(rigidbody.get_body(0).get_atom(0).coordinates() == Vector3<double>( 1, -3, -1));
        CHECK(rigidbody.get_body(0).get_atom(1).coordinates() == Vector3<double>(-1, -3, -1));
        transform.undo();
        CHECK(rigidbody.get_body(0).get_atom(0).coordinates() == Vector3<double>(-1, -1, -1));
        CHECK(rigidbody.get_body(0).get_atom(1).coordinates() == Vector3<double>(-1,  1, -1));

        // multi-body translate
        transform.apply({Matrix<double>::identity(3), Vector3<double>(1, 0, 0)}, manager->distance_constraints[1]);
        CHECK(rigidbody.get_body(0).get_atom(0).coordinates() == Vector3<double>(0, -1, -1));
        CHECK(rigidbody.get_body(0).get_atom(1).coordinates() == Vector3<double>(0,  1, -1));
        CHECK(rigidbody.get_body(1).get_atom(0).coordinates() == Vector3<double>(2, -1, -1));
        CHECK(rigidbody.get_body(1).get_atom(1).coordinates() == Vector3<double>(2,  1, -1));
        transform.undo();
        CHECK(rigidbody.get_body(0).get_atom(0).coordinates() == Vector3<double>(-1, -1, -1));
        CHECK(rigidbody.get_body(0).get_atom(1).coordinates() == Vector3<double>(-1,  1, -1));
        CHECK(rigidbody.get_body(1).get_atom(0).coordinates() == Vector3<double>( 1, -1, -1));
        CHECK(rigidbody.get_body(1).get_atom(1).coordinates() == Vector3<double>( 1,  1, -1));

        // multi-body rotate
        transform.apply({matrix::rotation_matrix(Vector3<double>{0, 0, 1}, M_PI/2), Vector3<double>(0, 0, 0)}, manager->distance_constraints[1]);
        CHECK(rigidbody.get_body(0).get_atom(0).coordinates() == Vector3<double>( 5, -3, -1));
        CHECK(rigidbody.get_body(0).get_atom(1).coordinates() == Vector3<double>( 3, -3, -1));
        CHECK(rigidbody.get_body(1).get_atom(0).coordinates() == Vector3<double>( 5, -1, -1));
        CHECK(rigidbody.get_body(1).get_atom(1).coordinates() == Vector3<double>( 3, -1, -1));
        transform.undo();
        CHECK(rigidbody.get_body(0).get_atom(0).coordinates() == Vector3<double>(-1, -1, -1));
        CHECK(rigidbody.get_body(0).get_atom(1).coordinates() == Vector3<double>(-1,  1, -1));
        CHECK(rigidbody.get_body(1).get_atom(0).coordinates() == Vector3<double>( 1, -1, -1));
        CHECK(rigidbody.get_body(1).get_atom(1).coordinates() == Vector3<double>( 1,  1, -1));
    }
}

auto vector_contains = [] (std::vector<unsigned int> vec, std::vector<unsigned int> vals) {
    std::unordered_set<unsigned int> set(vec.begin(), vec.end());
    for (auto val : vals) {
        if (!set.contains(val)) {
            return false;
        }
    }
    return true;
};

TEST_CASE_METHOD(fixture, "RigidTransform::get_connected") {
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::None;
    settings::general::verbose = false;

    SECTION("get_connected") {
        struct TestRigidTransform : public transform::RigidTransform {
            using RigidTransform::RigidTransform;
            using RigidTransform::get_connected;
        };

        SECTION("simple") {
            RigidBody rigidbody(bodies);
            auto manager = rigidbody.get_constraint_manager();

            manager->add_constraint(rigidbody::constraints::DistanceConstraint(&rigidbody, 0, 1, 0, 0)); // 0
            manager->add_constraint(rigidbody::constraints::DistanceConstraint(&rigidbody, 1, 2, 0, 0)); // 1
            manager->add_constraint(rigidbody::constraints::DistanceConstraint(&rigidbody, 2, 3, 0, 0)); // 2
            manager->add_constraint(rigidbody::constraints::DistanceConstraint(&rigidbody, 3, 4, 0, 0)); // 3

            TestRigidTransform transform(&rigidbody);

            // 0 - 1 - 2 - 3 - 4            //
            auto group1 = transform.get_connected(rigidbody.get_constraint_manager()->distance_constraints[0]);
            REQUIRE(group1.indices.size() == 1);
            CHECK(group1.indices[0] == 0);

            auto group2 = transform.get_connected(rigidbody.get_constraint_manager()->distance_constraints[1]);
            REQUIRE(group2.indices.size() == 2);
            CHECK(vector_contains(group2.indices, {0, 1}));

            auto group3 = transform.get_connected(rigidbody.get_constraint_manager()->distance_constraints[2]);
            REQUIRE(group3.indices.size() == 2);
            CHECK(vector_contains(group3.indices, {3, 4}));

            auto group4 = transform.get_connected(rigidbody.get_constraint_manager()->distance_constraints[3]);
            REQUIRE(group4.indices.size() == 1);
            CHECK(group4.indices[0] == 4);
        }

        SECTION("complex") {
            Body b6(std::vector<AtomFF>{AtomFF({0, 1, 0}, form_factor::form_factor_t::C), AtomFF({0, 2, 0}, form_factor::form_factor_t::C)});
            Body b7(std::vector<AtomFF>{AtomFF({0, 3, 0}, form_factor::form_factor_t::C), AtomFF({0, 4, 0}, form_factor::form_factor_t::C)});
            Body b8(std::vector<AtomFF>{AtomFF({1, 0, 0}, form_factor::form_factor_t::C), AtomFF({2, 0, 0}, form_factor::form_factor_t::C)});
            bodies = {b1, b2, b3, b4, b5, b6, b7, b8};
            RigidBody rigidbody(bodies);
            auto manager = rigidbody.get_constraint_manager();

            manager->add_constraint(rigidbody::constraints::DistanceConstraint(&rigidbody, 0, 1, 0, 0)); // 0
            manager->add_constraint(rigidbody::constraints::DistanceConstraint(&rigidbody, 1, 2, 0, 0)); // 1
            manager->add_constraint(rigidbody::constraints::DistanceConstraint(&rigidbody, 2, 3, 0, 0)); // 2
            manager->add_constraint(rigidbody::constraints::DistanceConstraint(&rigidbody, 3, 4, 0, 0)); // 3
            manager->add_constraint(rigidbody::constraints::DistanceConstraint(&rigidbody, 3, 5, 0, 0)); // 4
            manager->add_constraint(rigidbody::constraints::DistanceConstraint(&rigidbody, 5, 6, 0, 0)); // 5
            manager->add_constraint(rigidbody::constraints::DistanceConstraint(&rigidbody, 3, 7, 0, 0)); // 6

            //             5 - 6            //
            //             |                //
            // 0 - 1 - 2 - 3 - 4            //
            //             |                //
            //             7                //

            TestRigidTransform transform(&rigidbody);
            auto group1 = transform.get_connected(rigidbody.get_constraint_manager()->distance_constraints[0]);
            REQUIRE(group1.indices.size() == 1);
            CHECK(group1.indices[0] == 0);

            auto group2 = transform.get_connected(rigidbody.get_constraint_manager()->distance_constraints[1]);
            REQUIRE(group2.indices.size() == 2);
            CHECK(vector_contains(group2.indices, {0, 1}));

            auto group3 = transform.get_connected(rigidbody.get_constraint_manager()->distance_constraints[2]);
            REQUIRE(group3.indices.size() == 3);
            CHECK(vector_contains(group3.indices, {0, 1, 2}));

            auto group4 = transform.get_connected(rigidbody.get_constraint_manager()->distance_constraints[3]);
            REQUIRE(group4.indices.size() == 1);
            CHECK(group4.indices[0] == 4);

            auto group5 = transform.get_connected(rigidbody.get_constraint_manager()->distance_constraints[4]);
            REQUIRE(group5.indices.size() == 2);
            CHECK(vector_contains(group5.indices, {5, 6}));

            auto group6 = transform.get_connected(rigidbody.get_constraint_manager()->distance_constraints[5]);
            REQUIRE(group6.indices.size() == 1);
            CHECK(group6.indices[0] == 6);

            auto group7 = transform.get_connected(rigidbody.get_constraint_manager()->distance_constraints[6]);
            REQUIRE(group7.indices.size() == 1);
            CHECK(group7.indices[0] == 7);
        }
    }
}