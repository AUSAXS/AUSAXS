#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <rigidbody/Rigidbody.h>
#include <rigidbody/transform/RigidTransform.h>
#include <rigidbody/constraints/DistanceConstraintBond.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/detail/SystemSpecification.h>
#include <data/Body.h>
#include <data/Molecule.h>
#include <math/MatrixUtils.h>
#include <settings/All.h>

#include <numbers>

using namespace ausaxs;
using namespace ausaxs::data;
using namespace ausaxs::rigidbody;

TEST_CASE("RigidTransform::apply single body group") {
    settings::molecule::implicit_hydrogens = false;
    settings::molecule::center = false;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::None;
    settings::grid::min_bins = 10;

    SECTION("single body behaves like SingleTransform") {
        AtomFF a1({0, 0, 0}, form_factor::form_factor_t::C);
        AtomFF a2({1, 0, 0}, form_factor::form_factor_t::C);
        
        Body b1(std::vector<AtomFF>{a1});
        Body b2(std::vector<AtomFF>{a2});
        
        Rigidbody rigidbody(Molecule{std::vector<Body>{b1, b2}});
        rigidbody.constraints->add_constraint(
            std::make_unique<constraints::DistanceConstraintBond>(&rigidbody.molecule, 0, 1)
        );
        
        transform::RigidTransform transformer(&rigidbody);
        auto constraint = rigidbody.constraints->discoverable_constraints[0].get();
        
        auto cm0_before = rigidbody.molecule.get_body(0).get_cm();
        auto cm1_before = rigidbody.molecule.get_body(1).get_cm();
        
        // Apply transformation
        transformer.apply({{0, 2, 0}, {0, 0, 0}}, constraint);
        
        // At least one body should have changed position
        auto cm0_after = rigidbody.molecule.get_body(0).get_cm();
        auto cm1_after = rigidbody.molecule.get_body(1).get_cm();
        
        bool body0_moved = (cm0_after - cm0_before).norm() > 1e-10;
        bool body1_moved = (cm1_after - cm1_before).norm() > 1e-10;
        REQUIRE((body0_moved || body1_moved));
    }
}

TEST_CASE("RigidTransform::apply multi-body group") {
    settings::molecule::implicit_hydrogens = false;
    settings::molecule::center = false;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::None;
    settings::grid::min_bins = 10;

    SECTION("linear chain - transform smaller side") {
        // Create chain: 0 - 1 - 2 - 3
        AtomFF a1({0, 0, 0}, form_factor::form_factor_t::C);
        AtomFF a2({1, 0, 0}, form_factor::form_factor_t::C);
        AtomFF a3({2, 0, 0}, form_factor::form_factor_t::C);
        AtomFF a4({3, 0, 0}, form_factor::form_factor_t::C);
        
        Body b1(std::vector<AtomFF>{a1});
        Body b2(std::vector<AtomFF>{a2});
        Body b3(std::vector<AtomFF>{a3});
        Body b4(std::vector<AtomFF>{a4});
        
        Rigidbody rigidbody(Molecule{std::vector<Body>{b1, b2, b3, b4}});
        rigidbody.constraints->add_constraint(
            std::make_unique<constraints::DistanceConstraintBond>(&rigidbody.molecule, 0, 1)
        );
        rigidbody.constraints->add_constraint(
            std::make_unique<constraints::DistanceConstraintBond>(&rigidbody.molecule, 1, 2)
        );
        rigidbody.constraints->add_constraint(
            std::make_unique<constraints::DistanceConstraintBond>(&rigidbody.molecule, 2, 3)
        );
        
        transform::RigidTransform transformer(&rigidbody);
        
        // Transform at constraint 1 (between bodies 1 and 2)
        // Should move the smaller group
        auto constraint = rigidbody.constraints->discoverable_constraints[1].get();
        
        auto cm0_before = rigidbody.molecule.get_body(0).get_cm();
        auto cm1_before = rigidbody.molecule.get_body(1).get_cm();
        auto cm2_before = rigidbody.molecule.get_body(2).get_cm();
        auto cm3_before = rigidbody.molecule.get_body(3).get_cm();
        
        transformer.apply({{0, 1, 0}, {0, 0, 0}}, constraint);
        
        // At least some bodies should have moved
        auto cm0_after = rigidbody.molecule.get_body(0).get_cm();
        auto cm1_after = rigidbody.molecule.get_body(1).get_cm();
        auto cm2_after = rigidbody.molecule.get_body(2).get_cm();
        auto cm3_after = rigidbody.molecule.get_body(3).get_cm();
        
        bool any_moved = (cm0_after - cm0_before).norm() > 1e-10 ||
                        (cm1_after - cm1_before).norm() > 1e-10 ||
                        (cm2_after - cm2_before).norm() > 1e-10 ||
                        (cm3_after - cm3_before).norm() > 1e-10;
        REQUIRE(any_moved);
    }

    SECTION("rigid rotation of group") {
        // Create simple two-body system
        AtomFF a1({-1, 0, 0}, form_factor::form_factor_t::C);
        AtomFF a2({0, 0, 0}, form_factor::form_factor_t::C);
        AtomFF a3({1, 0, 0}, form_factor::form_factor_t::C);
        
        Body b1(std::vector<AtomFF>{a1});
        Body b2(std::vector<AtomFF>{a2});
        Body b3(std::vector<AtomFF>{a3});
        
        Rigidbody rigidbody(Molecule{std::vector<Body>{b1, b2, b3}});
        rigidbody.constraints->add_constraint(
            std::make_unique<constraints::DistanceConstraintBond>(&rigidbody.molecule, 0, 1)
        );
        rigidbody.constraints->add_constraint(
            std::make_unique<constraints::DistanceConstraintBond>(&rigidbody.molecule, 1, 2)
        );
        
        transform::RigidTransform transformer(&rigidbody);
        
        // Transform at constraint 1 - should rotate bodies 0 and 1 together
        auto constraint = rigidbody.constraints->discoverable_constraints[1].get();
        
        // Store initial distance between bodies 0 and 1
        auto initial_dist = (rigidbody.molecule.get_body(0).get_cm() - 
                            rigidbody.molecule.get_body(1).get_cm()).norm();
        
        // Apply 90-degree rotation
        transformer.apply({{0, 0, 0}, {0, 0, std::numbers::pi/2}}, constraint);
        
        // Distance between bodies 0 and 1 should be preserved (rigid)
        auto final_dist = (rigidbody.molecule.get_body(0).get_cm() - 
                          rigidbody.molecule.get_body(1).get_cm()).norm();
        
        REQUIRE_THAT(final_dist, Catch::Matchers::WithinAbs(initial_dist, 1e-10));
    }
}

TEST_CASE("RigidTransform::apply branched structure") {
    settings::molecule::implicit_hydrogens = false;
    settings::molecule::center = false;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::None;
    settings::grid::scaling = 2;

    SECTION("T-shaped structure") {
        // Create T-shape:
        //     2
        //     |
        // 0 - 1 - 3
        
        AtomFF a1({-1, 0, 0}, form_factor::form_factor_t::C);
        AtomFF a2({0, 0, 0}, form_factor::form_factor_t::C);
        AtomFF a3({0, 1, 0}, form_factor::form_factor_t::C);
        AtomFF a4({1, 0, 0}, form_factor::form_factor_t::C);
        
        Body b1(std::vector<AtomFF>{a1});
        Body b2(std::vector<AtomFF>{a2});
        Body b3(std::vector<AtomFF>{a3});
        Body b4(std::vector<AtomFF>{a4});
        
        Rigidbody rigidbody(Molecule{std::vector<Body>{b1, b2, b3, b4}});
        rigidbody.constraints->add_constraint(
            std::make_unique<constraints::DistanceConstraintBond>(&rigidbody.molecule, 0, 1)
        );
        rigidbody.constraints->add_constraint(
            std::make_unique<constraints::DistanceConstraintBond>(&rigidbody.molecule, 1, 2)
        );
        rigidbody.constraints->add_constraint(
            std::make_unique<constraints::DistanceConstraintBond>(&rigidbody.molecule, 1, 3)
        );
        
        transform::RigidTransform transformer(&rigidbody);
        
        // Transform at constraint 0 - should move only body 0
        auto constraint0 = rigidbody.constraints->discoverable_constraints[0].get();
        transformer.apply({{0, 2, 0}, {0, 0, 0}}, constraint0);
        
        // Body 0 moved
        REQUIRE_THAT(rigidbody.molecule.get_body(0).get_cm().y(), Catch::Matchers::WithinAbs(2.0, 1e-10));
        
        // Others should not have moved
        REQUIRE_THAT(rigidbody.molecule.get_body(1).get_cm().x(), Catch::Matchers::WithinAbs(0.0, 1e-10));
        REQUIRE_THAT(rigidbody.molecule.get_body(2).get_cm().y(), Catch::Matchers::WithinAbs(1.0, 1e-10));
        REQUIRE_THAT(rigidbody.molecule.get_body(3).get_cm().x(), Catch::Matchers::WithinAbs(1.0, 1e-10));
    }
}

TEST_CASE("RigidTransform::undo") {
    settings::molecule::implicit_hydrogens = false;
    settings::molecule::center = false;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::None;
    settings::grid::min_bins = 25;

    SECTION("undo restores all bodies in group") {
        AtomFF a1({0, 0, 0}, form_factor::form_factor_t::C);
        AtomFF a2({1, 0, 0}, form_factor::form_factor_t::C);
        AtomFF a3({2, 0, 0}, form_factor::form_factor_t::C);
        
        Body b1(std::vector<AtomFF>{a1});
        Body b2(std::vector<AtomFF>{a2});
        Body b3(std::vector<AtomFF>{a3});
        
        Rigidbody rigidbody(Molecule{std::vector<Body>{b1, b2, b3}});
        rigidbody.constraints->add_constraint(
            std::make_unique<constraints::DistanceConstraintBond>(&rigidbody.molecule, 0, 1)
        );
        rigidbody.constraints->add_constraint(
            std::make_unique<constraints::DistanceConstraintBond>(&rigidbody.molecule, 1, 2)
        );
        
        transform::RigidTransform transformer(&rigidbody);
        auto constraint = rigidbody.constraints->discoverable_constraints[1].get();
        
        auto cm0_before = rigidbody.molecule.get_body(0).get_cm();
        auto cm1_before = rigidbody.molecule.get_body(1).get_cm();
        auto params0_before = rigidbody.conformation->absolute_parameters.parameters[0];
        auto params1_before = rigidbody.conformation->absolute_parameters.parameters[1];
        
        // Apply transformation
        transformer.apply({{5, 5, 5}, {0.5, 0.5, 0.5}}, constraint);
        
        // Undo
        transformer.undo();
        
        // Check restoration
        auto cm0_after = rigidbody.molecule.get_body(0).get_cm();
        auto cm1_after = rigidbody.molecule.get_body(1).get_cm();
        auto params0_after = rigidbody.conformation->absolute_parameters.parameters[0];
        auto params1_after = rigidbody.conformation->absolute_parameters.parameters[1];
        
        REQUIRE_THAT(cm0_after.x(), Catch::Matchers::WithinAbs(cm0_before.x(), 1e-10));
        REQUIRE_THAT(cm1_after.x(), Catch::Matchers::WithinAbs(cm1_before.x(), 1e-10));
        REQUIRE_THAT(params0_after.translation.x(), Catch::Matchers::WithinAbs(params0_before.translation.x(), 1e-10));
        REQUIRE_THAT(params1_after.translation.x(), Catch::Matchers::WithinAbs(params1_before.translation.x(), 1e-10));
    }
}
