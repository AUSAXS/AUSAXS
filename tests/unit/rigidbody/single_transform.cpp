#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <rigidbody/Rigidbody.h>
#include <rigidbody/transform/SingleTransform.h>
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

TEST_CASE("SingleTransform::apply basic transformations") {
    settings::molecule::implicit_hydrogens = false;
    settings::molecule::center = false;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::None;
    settings::grid::scaling = 2;

    SECTION("translation only") {
        AtomFF a1({0, 0, 0}, form_factor::form_factor_t::C);
        AtomFF a2({1, 0, 0}, form_factor::form_factor_t::C);
        
        Body b1(std::vector<AtomFF>{a1});
        Body b2(std::vector<AtomFF>{a2});
        
        Rigidbody rigidbody(Molecule{std::vector<Body>{b1, b2}});
        rigidbody.constraints->add_constraint(
            std::make_unique<constraints::DistanceConstraintBond>(&rigidbody.molecule, 0, 1)
        );
        
        transform::SingleTransform transformer(&rigidbody);
        auto constraint = rigidbody.constraints->discoverable_constraints[0].get();
        
        // Apply translation
        transformer.apply({{0, 1, 0}, {0, 0, 0}}, constraint);
        
        // Only body 0 should have moved
        REQUIRE_THAT(rigidbody.molecule.get_body(0).get_cm().y(), Catch::Matchers::WithinAbs(1.0, 1e-10));
        // Body 1 should not have moved
        REQUIRE_THAT(rigidbody.molecule.get_body(1).get_cm().x(), Catch::Matchers::WithinAbs(1.0, 1e-10));
    }

    SECTION("rotation only") {
        AtomFF a1({-1, 0, 0}, form_factor::form_factor_t::C);
        AtomFF a2({0, 0, 0}, form_factor::form_factor_t::C);
        
        Body b1(std::vector<AtomFF>{a1});
        Body b2(std::vector<AtomFF>{a2});
        
        Rigidbody rigidbody(Molecule{std::vector<Body>{b1, b2}});
        rigidbody.constraints->add_constraint(
            std::make_unique<constraints::DistanceConstraintBond>(&rigidbody.molecule, 0, 1)
        );
        
        transform::SingleTransform transformer(&rigidbody);
        auto constraint = rigidbody.constraints->discoverable_constraints[0].get();
        
        // Apply 90-degree rotation around z-axis
        transformer.apply({{0, 0, 0}, {0, 0, std::numbers::pi/2}}, constraint);
        
        // Body should have rotated around the pivot (body 1's position)
        auto cm = rigidbody.molecule.get_body(0).get_cm();
        REQUIRE_THAT(cm.x(), Catch::Matchers::WithinAbs(0.0, 1e-10));
        REQUIRE_THAT(cm.y(), Catch::Matchers::WithinAbs(-1.0, 1e-10));
    }

    SECTION("combined rotation and translation") {
        AtomFF a1({1, 0, 0}, form_factor::form_factor_t::C);
        AtomFF a2({0, 0, 0}, form_factor::form_factor_t::C);
        
        Body b1(std::vector<AtomFF>{a1});
        Body b2(std::vector<AtomFF>{a2});
        
        Rigidbody rigidbody(Molecule{std::vector<Body>{b1, b2}});
        rigidbody.constraints->add_constraint(
            std::make_unique<constraints::DistanceConstraintBond>(&rigidbody.molecule, 0, 1)
        );
        
        transform::SingleTransform transformer(&rigidbody);
        auto constraint = rigidbody.constraints->discoverable_constraints[0].get();
        
        // Apply rotation then translation
        transformer.apply({{1, 1, 0}, {0, 0, std::numbers::pi/2}}, constraint);
        
        // Body rotated 90 degrees then translated
        auto cm = rigidbody.molecule.get_body(0).get_cm();
        REQUIRE_THAT(cm.x(), Catch::Matchers::WithinAbs(1.0, 1e-10));
        REQUIRE_THAT(cm.y(), Catch::Matchers::WithinAbs(2.0, 1e-10));
    }
}

TEST_CASE("SingleTransform::undo") {
    settings::molecule::implicit_hydrogens = false;
    settings::molecule::center = false;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::None;
    settings::grid::scaling = 2;

    SECTION("undo restores original state") {
        AtomFF a1({1, 1, 1}, form_factor::form_factor_t::C);
        AtomFF a2({0, 0, 0}, form_factor::form_factor_t::C);
        
        Body b1(std::vector<AtomFF>{a1});
        Body b2(std::vector<AtomFF>{a2});
        
        Rigidbody rigidbody(Molecule{std::vector<Body>{b1, b2}});
        rigidbody.constraints->add_constraint(
            std::make_unique<constraints::DistanceConstraintBond>(&rigidbody.molecule, 0, 1)
        );
        
        transform::SingleTransform transformer(&rigidbody);
        auto constraint = rigidbody.constraints->discoverable_constraints[0].get();
        
        auto original_cm = rigidbody.molecule.get_body(0).get_cm();
        auto original_params = rigidbody.conformation->absolute_parameters.parameters[0];
        
        // Apply transformation
        transformer.apply({{2, 3, 4}, {0.1, 0.2, 0.3}}, constraint);
        
        // Verify it changed
        auto changed_cm = rigidbody.molecule.get_body(0).get_cm();
        REQUIRE(changed_cm != original_cm);
        
        // Undo
        transformer.undo();
        
        // Verify restoration
        auto restored_cm = rigidbody.molecule.get_body(0).get_cm();
        REQUIRE_THAT(restored_cm.x(), Catch::Matchers::WithinAbs(original_cm.x(), 1e-10));
        REQUIRE_THAT(restored_cm.y(), Catch::Matchers::WithinAbs(original_cm.y(), 1e-10));
        REQUIRE_THAT(restored_cm.z(), Catch::Matchers::WithinAbs(original_cm.z(), 1e-10));
        
        auto& restored_params = rigidbody.conformation->absolute_parameters.parameters[0];
        REQUIRE_THAT(restored_params.translation.x(), Catch::Matchers::WithinAbs(original_params.translation.x(), 1e-10));
        REQUIRE_THAT(restored_params.rotation.x(), Catch::Matchers::WithinAbs(original_params.rotation.x(), 1e-10));
    }
}

TEST_CASE("SingleTransform::apply parameter reconstruction") {
    settings::molecule::implicit_hydrogens = false;
    settings::molecule::center = false;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::None;
    settings::grid::scaling = 2;

    SECTION("parameters can reconstruct body state") {
        AtomFF a1({3, 2, 1}, form_factor::form_factor_t::C);
        AtomFF a2({0, 0, 0}, form_factor::form_factor_t::C);
        
        Body b1(std::vector<AtomFF>{a1});
        Body b2(std::vector<AtomFF>{a2});
        
        Rigidbody rigidbody(Molecule{std::vector<Body>{b1, b2}});
        rigidbody.constraints->add_constraint(
            std::make_unique<constraints::DistanceConstraintBond>(&rigidbody.molecule, 0, 1)
        );
        
        transform::SingleTransform transformer(&rigidbody);
        auto constraint = rigidbody.constraints->discoverable_constraints[0].get();
        
        // Apply transformation
        transformer.apply({{1, 2, 3}, {0.5, 0.3, 0.1}}, constraint);
        
        // Reconstruct from parameters
        auto& current_body = rigidbody.molecule.get_body(0);
        auto& params = rigidbody.conformation->absolute_parameters.parameters[0];
        auto& original = rigidbody.conformation->initial_conformation[0];
        
        Body reconstructed = original;
        reconstructed.rotate(matrix::rotation_matrix(params.rotation));
        reconstructed.translate(params.translation);
        
        auto current_cm = current_body.get_cm();
        auto reconstructed_cm = reconstructed.get_cm();
        
        REQUIRE_THAT(reconstructed_cm.x(), Catch::Matchers::WithinAbs(current_cm.x(), 1e-3));
        REQUIRE_THAT(reconstructed_cm.y(), Catch::Matchers::WithinAbs(current_cm.y(), 1e-3));
        REQUIRE_THAT(reconstructed_cm.z(), Catch::Matchers::WithinAbs(current_cm.z(), 1e-3));
    }
}

TEST_CASE("SingleTransform::apply multiple sequential transformations") {
    settings::molecule::implicit_hydrogens = false;
    settings::molecule::center = false;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::None;
    settings::grid::scaling = 2;

    SECTION("multiple transformations accumulate correctly") {
        AtomFF a1({1, 0, 0}, form_factor::form_factor_t::C);
        AtomFF a2({0, 0, 0}, form_factor::form_factor_t::C);
        
        Body b1(std::vector<AtomFF>{a1});
        Body b2(std::vector<AtomFF>{a2});
        
        Rigidbody rigidbody(Molecule{std::vector<Body>{b1, b2}});
        rigidbody.constraints->add_constraint(
            std::make_unique<constraints::DistanceConstraintBond>(&rigidbody.molecule, 0, 1)
        );
        
        transform::SingleTransform transformer(&rigidbody);
        auto constraint = rigidbody.constraints->discoverable_constraints[0].get();
        
        // Apply first transformation
        transformer.apply({{0, 0, 0}, {0, 0, std::numbers::pi/4}}, constraint);
        
        // Apply second transformation
        transformer.apply({{0, 0, 0}, {0, 0, std::numbers::pi/4}}, constraint);
        
        // Total rotation should be pi/2
        auto& params = rigidbody.conformation->absolute_parameters.parameters[0];
        REQUIRE_THAT(params.rotation.z(), Catch::Matchers::WithinAbs(std::numbers::pi/2, 1e-6));
        
        // Position should be correct
        auto cm = rigidbody.molecule.get_body(0).get_cm();
        REQUIRE_THAT(cm.x(), Catch::Matchers::WithinAbs(0.0, 1e-10));
        REQUIRE_THAT(cm.y(), Catch::Matchers::WithinAbs(1.0, 1e-10));
    }
}

TEST_CASE("SingleTransform::apply only affects single body") {
    settings::molecule::implicit_hydrogens = false;
    settings::molecule::center = false;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::None;
    settings::grid::scaling = 2;

    SECTION("other bodies remain unchanged") {
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
        
        transform::SingleTransform transformer(&rigidbody);
        auto constraint = rigidbody.constraints->discoverable_constraints[0].get();
        
        auto body1_cm_before = rigidbody.molecule.get_body(1).get_cm();
        auto body2_cm_before = rigidbody.molecule.get_body(2).get_cm();
        
        // Transform body 0
        transformer.apply({{5, 5, 5}, {0, 0, 0}}, constraint);
        
        // Bodies 1 and 2 should not have moved
        auto body1_cm_after = rigidbody.molecule.get_body(1).get_cm();
        auto body2_cm_after = rigidbody.molecule.get_body(2).get_cm();
        
        REQUIRE_THAT(body1_cm_after.x(), Catch::Matchers::WithinAbs(body1_cm_before.x(), 1e-10));
        REQUIRE_THAT(body1_cm_after.y(), Catch::Matchers::WithinAbs(body1_cm_before.y(), 1e-10));
        REQUIRE_THAT(body2_cm_after.x(), Catch::Matchers::WithinAbs(body2_cm_before.x(), 1e-10));
        REQUIRE_THAT(body2_cm_after.y(), Catch::Matchers::WithinAbs(body2_cm_before.y(), 1e-10));
    }
}
