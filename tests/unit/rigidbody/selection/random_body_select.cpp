#include <catch2/catch_test_macros.hpp>

#include <data/Body.h>
#include <data/Molecule.h>
#include <rigidbody/Rigidbody.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/constraints/DistanceConstraintCM.h>
#include <rigidbody/selection/RandomBodySelect.h>
#include <settings/All.h>

using namespace ausaxs;
using namespace ausaxs::data;
using namespace ausaxs::rigidbody;
using namespace ausaxs::rigidbody::selection;

struct SelectionStrategiesFixture {
    SelectionStrategiesFixture() {
        settings::general::verbose = false;
        settings::molecule::implicit_hydrogens = false;
        
        // Create test bodies
        AtomFF a1({0, 0, 0}, form_factor::form_factor_t::C);
        AtomFF a2({5, 0, 0}, form_factor::form_factor_t::C);
        AtomFF a3({10, 0, 0}, form_factor::form_factor_t::C);
        AtomFF a4({15, 0, 0}, form_factor::form_factor_t::C);
        
        Body b1(std::vector{a1});
        Body b2(std::vector{a2});
        Body b3(std::vector{a3});
        Body b4(std::vector{a4});
        
        rb = std::make_unique<Rigidbody>(Molecule{std::vector<Body>{b1, b2, b3, b4}});
        
        // Add some constraints for constraint-based selection
        rb->constraints->add_constraint(
            std::make_unique<constraints::DistanceConstraintCM>(&rb->molecule, 0, 1)
        );
        rb->constraints->add_constraint(
            std::make_unique<constraints::DistanceConstraintCM>(&rb->molecule, 1, 2)
        );
        rb->constraints->add_constraint(
            std::make_unique<constraints::DistanceConstraintCM>(&rb->molecule, 2, 3)
        );
    }
    
    std::unique_ptr<Rigidbody> rb;
};

TEST_CASE_METHOD(SelectionStrategiesFixture, "SelectionStrategies::RandomBodySelect") {
    RandomBodySelect selector(rb.get());
    
    for (int i = 0; i < 10; ++i) {
        auto [ibody, iconstraint, isymmetry] = selector.next(ParameterMask::all());
        CHECK(ibody < static_cast<int>(rb->molecule.get_bodies().size()));
    }
}

// Regression test: strategies used to cache the body count at construction. If bodies were later removed (e.g. by the "delete", "merge", or "convert_to_symmetry" 
// sequencer elements), next() could still return indices from the old, larger range, and callers like LimitedParameterGenerator::next() would index out of bounds.
TEST_CASE_METHOD(SelectionStrategiesFixture, "RandomBodySelect: stale body count after removal") {
    RandomBodySelect selector(rb.get());
    rb->molecule.get_bodies().resize(1);

    for (int i = 0; i < 50; ++i) {
        auto [ibody, iconstraint, isymmetry] = selector.next(ParameterMask::all());
        CHECK(ibody < static_cast<int>(rb->molecule.get_bodies().size()));
    }
}
