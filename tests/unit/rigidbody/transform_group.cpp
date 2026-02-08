#include <catch2/catch_test_macros.hpp>

#include <rigidbody/transform/TransformGroup.h>
#include <rigidbody/constraints/DistanceConstraint.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <data/atoms/AtomFF.h>
#include <math/Vector3.h>
#include <settings/All.h>
#include <form_factor/FormFactorType.h>

using namespace ausaxs;
using namespace ausaxs::rigidbody;
using namespace ausaxs::rigidbody::transform;
using namespace ausaxs::rigidbody::constraints;
using namespace ausaxs::data;

TEST_CASE("TransformGroup::construction") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::molecule::center = false;
    
    AtomFF a1 = AtomFF({-1, -1, -1}, form_factor::form_factor_t::C);
    AtomFF a2 = AtomFF({-1,  1, -1}, form_factor::form_factor_t::C);
    AtomFF a3 = AtomFF({-1, -1,  1}, form_factor::form_factor_t::C);
    AtomFF a4 = AtomFF({-1,  1,  1}, form_factor::form_factor_t::C);
    
    Body b1 = Body(std::vector{a1, a2});
    Body b2 = Body(std::vector{a3, a4});
    std::vector<Body> bodies_vec = {b1, b2};
    Molecule molecule(bodies_vec);
    auto& bodies = molecule.get_bodies();
    
    SECTION("construct with single body") {
        DistanceConstraint constraint(&molecule, 0, 1, 0, 0);
        std::vector<data::Body*> body_ptrs = {&bodies[0]};
        std::vector<unsigned int> indices = {0};
        Vector3<double> pivot(0, 0, 0);
        
        TransformGroup group(body_ptrs, indices, constraint, pivot);
        
        CHECK(group.bodies.size() == 1);
        CHECK(group.indices.size() == 1);
        CHECK(group.indices[0] == 0);
        CHECK(group.pivot.x() == 0);
        CHECK(group.pivot.y() == 0);
        CHECK(group.pivot.z() == 0);
    }

    SECTION("construct with multiple bodies") {
        DistanceConstraint constraint(&molecule, 0, 1, 0, 0);
        std::vector<data::Body*> body_ptrs = {&bodies[0], &bodies[1]};
        std::vector<unsigned int> indices = {0, 1};
        Vector3<double> pivot(1.5, 2.5, 3.5);
        
        TransformGroup group(body_ptrs, indices, constraint, pivot);
        
        CHECK(group.bodies.size() == 2);
        CHECK(group.indices.size() == 2);
        CHECK(group.indices[0] == 0);
        CHECK(group.indices[1] == 1);
        CHECK(group.pivot.x() == 1.5);
        CHECK(group.pivot.y() == 2.5);
        CHECK(group.pivot.z() == 3.5);
    }
}
