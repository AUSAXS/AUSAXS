#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <rigidbody/transform/TransformGroup.h>
#include <rigidbody/constraints/DistanceConstraintBond.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/Rigidbody.h>
#include <data/Body.h>
#include <data/Molecule.h>
#include <settings/All.h>

using namespace ausaxs;
using namespace ausaxs::data;
using namespace ausaxs::rigidbody;

TEST_CASE("TransformGroup::TransformGroup") {
    settings::molecule::implicit_hydrogens = false;
    settings::molecule::center = false;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::None;
    settings::grid::scaling = 2;

    SECTION("construction with single body") {
        AtomFF a1({0, 0, 0}, form_factor::form_factor_t::C);
        Body b1(std::vector<AtomFF>{a1});
        
        Rigidbody rigidbody(Molecule{std::vector<Body>{b1}});
        
        std::vector<observer_ptr<data::Body>> bodies = {&rigidbody.molecule.get_body(0)};
        std::vector<unsigned int> indices = {0};
        Vector3<double> pivot(1, 2, 3);
        
        transform::TransformGroup group(bodies, indices, nullptr, pivot);
        
        CHECK(group.bodies.size() == 1);
        CHECK(group.indices.size() == 1);
        CHECK(group.indices[0] == 0);
        CHECK(group.pivot.x() == 1);
        CHECK(group.pivot.y() == 2);
        CHECK(group.pivot.z() == 3);
    }

    SECTION("construction with multiple bodies") {
        AtomFF a1({0, 0, 0}, form_factor::form_factor_t::C);
        AtomFF a2({1, 0, 0}, form_factor::form_factor_t::C);
        AtomFF a3({2, 0, 0}, form_factor::form_factor_t::C);
        
        Body b1(std::vector<AtomFF>{a1});
        Body b2(std::vector<AtomFF>{a2});
        Body b3(std::vector<AtomFF>{a3});
        
        Rigidbody rigidbody(Molecule{std::vector<Body>{b1, b2, b3}});
        
        std::vector<observer_ptr<data::Body>> bodies = {
            &rigidbody.molecule.get_body(0),
            &rigidbody.molecule.get_body(1),
            &rigidbody.molecule.get_body(2)
        };
        std::vector<unsigned int> indices = {0, 1, 2};
        Vector3<double> pivot(0, 0, 0);
        
        transform::TransformGroup group(bodies, indices, nullptr, pivot);
        
        CHECK(group.bodies.size() == 3);
        CHECK(group.indices.size() == 3);
        CHECK(group.indices[0] == 0);
        CHECK(group.indices[1] == 1);
        CHECK(group.indices[2] == 2);
    }

    SECTION("construction with constraint") {
        settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::None;
        
        AtomFF a1({0, 0, 0}, form_factor::form_factor_t::C);
        AtomFF a2({1, 0, 0}, form_factor::form_factor_t::C);
        
        Body b1(std::vector<AtomFF>{a1});
        Body b2(std::vector<AtomFF>{a2});
        
        Rigidbody rigidbody(Molecule{std::vector<Body>{b1, b2}});
        rigidbody.constraints->add_constraint(
            std::make_unique<constraints::DistanceConstraintBond>(&rigidbody.molecule, 0, 1)
        );
        
        auto constraint = rigidbody.constraints->discoverable_constraints[0].get();
        std::vector<observer_ptr<data::Body>> bodies = {&rigidbody.molecule.get_body(0)};
        std::vector<unsigned int> indices = {0};
        Vector3<double> pivot = constraint->get_atom1().coordinates();
        
        transform::TransformGroup group(bodies, indices, constraint, pivot);
        
        CHECK(group.target == constraint);
        CHECK(group.pivot == constraint->get_atom1().coordinates());
    }
}
