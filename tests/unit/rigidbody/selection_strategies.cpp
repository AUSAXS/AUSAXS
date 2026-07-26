#include <catch2/catch_test_macros.hpp>

#include <rigidbody/selection/RandomBodySelect.h>
#include <rigidbody/selection/SequentialBodySelect.h>
#include <rigidbody/selection/RandomConstraintSelect.h>
#include <rigidbody/selection/SequentialConstraintSelect.h>
#include <rigidbody/selection/ManualSelect.h>
#include <rigidbody/Rigidbody.h>
#include <rigidbody/BodySplitter.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/constraints/DistanceConstraintCM.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <settings/All.h>

using namespace ausaxs;
using namespace ausaxs::data;
using namespace ausaxs::rigidbody;
using namespace ausaxs::rigidbody::selection;

struct SelectionStrategiesFixture {
    SelectionStrategiesFixture() {
        settings::general::verbose = false;
        settings::molecule::implicit_hydrogens = false;
        settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::None;
        
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
    
    SECTION("next returns valid body indices") {
        for (int i = 0; i < 10; ++i) {
            auto [ibody, iconstraint] = selector.next();
            CHECK(ibody < rb->molecule.get_bodies().size());
        }
    }
}

TEST_CASE_METHOD(SelectionStrategiesFixture, "SelectionStrategies::ManualSelect") {
    SECTION("next always returns the configured body") {
        ManualSelect selector(rb.get(), 2);

        for (int i = 0; i < 10; ++i) {
            auto [ibody, iconstraint] = selector.next();
            CHECK(ibody == 2);
        }
    }

    SECTION("next returns -1 for a body with no constraints") {
        Rigidbody isolated(Molecule{std::vector<Body>{
            Body(std::vector{AtomFF({0, 0, 0}, form_factor::form_factor_t::C)}),
            Body(std::vector{AtomFF({5, 0, 0}, form_factor::form_factor_t::C)})
        }});
        ManualSelect selector(&isolated, 0);

        auto [ibody, iconstraint] = selector.next();
        CHECK(ibody == 0);
        CHECK(iconstraint == -1);
    }
}

TEST_CASE_METHOD(SelectionStrategiesFixture, "SelectionStrategies::SequentialBodySelect") {
    SequentialBodySelect selector(rb.get());

    SECTION("next cycles through bodies in order") {
        unsigned int num_bodies = rb->molecule.get_bodies().size();

        for (unsigned int i = 0; i < num_bodies * 2; ++i) {
            auto [ibody, iconstraint] = selector.next();
            CHECK(ibody == i % num_bodies);
        }
    }
}

// Regression test: strategies used to cache the body count at construction. If bodies were later removed (e.g. by the "delete", "merge", or "convert_to_symmetry" 
// sequencer elements), next() could still return indices from the old, larger range, and callers like LimitedParameterGenerator::next() would index out of bounds.
TEST_CASE_METHOD(SelectionStrategiesFixture, "SelectionStrategies::stale body count after removal") {
    SECTION("RandomBodySelect") {
        RandomBodySelect selector(rb.get());
        rb->molecule.get_bodies().resize(1);

        for (int i = 0; i < 50; ++i) {
            auto [ibody, iconstraint] = selector.next();
            CHECK(ibody < rb->molecule.get_bodies().size());
        }
    }

    SECTION("SequentialBodySelect") {
        SequentialBodySelect selector(rb.get());
        rb->molecule.get_bodies().resize(1);

        for (int i = 0; i < 50; ++i) {
            auto [ibody, iconstraint] = selector.next();
            CHECK(ibody < rb->molecule.get_bodies().size());
        }
    }

    SECTION("RandomConstraintSelect") {
        // DistanceConstraintBond (the discoverable constraint type RandomConstraintSelect samples from) requires C-alpha backbone metadata, so load a real 
        // structure via BodySplitter with the Backbone generation strategy instead of hand-building one (see distance_constraint_bond.cpp)
        settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::Backbone;
        Rigidbody local_rb = BodySplitter::split("tests/files/LAR1-4.pdb", {9, 99, 202, 292});
        REQUIRE(local_rb.constraints->discoverable_constraints.size() == 4); // one bond between each of the 5 sequential bodies

        RandomConstraintSelect selector(&local_rb);

        // shrink the constraint list after the selector is constructed; next() used to sample from a distribution fixed at construction time, so this could 
        // pick a since-removed constraint
        local_rb.constraints->discoverable_constraints.resize(1);

        for (int i = 0; i < 50; ++i) {
            auto [ibody, iconstraint] = selector.next();
            CHECK(ibody < local_rb.molecule.get_bodies().size());
        }

        settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::None;
    }

    SECTION("RandomConstraintSelect throws instead of crashing once no constraints remain") {
        RandomConstraintSelect selector(rb.get());

        // the fixture never registers any discoverable constraints (only non-discoverable DistanceConstraintCM ones), so next() must throw rather than build 
        // a distribution over an empty range (UB from uniform_int_distribution(0, -1))
        REQUIRE(rb->constraints->discoverable_constraints.empty());
        CHECK_THROWS(selector.next());
    }

    SECTION("SequentialConstraintSelect") {
        SequentialConstraintSelect selector(rb.get());
        rb->molecule.get_bodies().resize(1);

        for (int i = 0; i < 50; ++i) {
            auto [ibody, iconstraint] = selector.next();
            CHECK(ibody < rb->molecule.get_bodies().size());
        }
    }
}
