#include <catch2/catch_test_macros.hpp>

#include <data/Body.h>
#include <data/Molecule.h>
#include <rigidbody/BodySplitter.h>
#include <rigidbody/Rigidbody.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/constraints/DistanceConstraintCM.h>
#include <rigidbody/selection/RandomConstraintSelect.h>
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

// Regression test: strategies used to cache the body count at construction. If bodies were later removed (e.g. by the "delete", "merge", or "convert_to_symmetry" 
// sequencer elements), next() could still return indices from the old, larger range, and callers like LimitedParameterGenerator::next() would index out of bounds.
TEST_CASE_METHOD(SelectionStrategiesFixture, "RandomConstraintSelect: stale body count after removal") {

    SECTION("RandomConstraintSelect") {
        // DistanceConstraintBond (the discoverable constraint type RandomConstraintSelect samples from) requires C-alpha backbone metadata, so load a real 
        // structure via BodySplitter and generate backbone constraints from it, instead of hand-building one (see distance_constraint_bond.cpp)
        Rigidbody local_rb = BodySplitter::split("tests/files/LAR1-4.pdb", {9, 99, 202, 292});
        local_rb.constraints->generate_constraints(settings::rigidbody::ConstraintGenerationStrategyChoice::Backbone);
        REQUIRE(local_rb.constraints->discoverable_constraints.size() == 4); // one bond between each of the 5 sequential bodies

        RandomConstraintSelect selector(&local_rb);

        // shrink the constraint list after the selector is constructed; next() used to sample from a distribution fixed at construction time, so this could 
        // pick a since-removed constraint
        local_rb.constraints->discoverable_constraints.resize(1);

        for (int i = 0; i < 50; ++i) {
            auto [ibody, iconstraint, isymmetry] = selector.next(ParameterMask::all());
            CHECK(ibody < static_cast<int>(local_rb.molecule.get_bodies().size()));
        }

    }

    SECTION("RandomConstraintSelect throws instead of crashing once no constraints remain") {
        RandomConstraintSelect selector(rb.get());

        // the fixture never registers any discoverable constraints (only non-discoverable DistanceConstraintCM ones), so next() must throw rather than build 
        // a distribution over an empty range (UB from uniform_int_distribution(0, -1))
        REQUIRE(rb->constraints->discoverable_constraints.empty());
        CHECK_THROWS(selector.next(ParameterMask::all()));
    }
}
