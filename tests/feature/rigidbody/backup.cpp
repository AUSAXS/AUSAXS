#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <rigidbody/Rigidbody.h>
#include <rigidbody/BodySplitter.h>
#include <rigidbody/detail/Conformation.h>
#include <rigidbody/detail/Configuration.h>
#include <rigidbody/parameters/BodyTransformParameters.h>
#include <rigidbody/parameters/RelativeTransformParameters.h>
#include <rigidbody/parameters/ParameterGenerationStrategy.h>
#include <rigidbody/transform/TransformStrategy.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/constraints/DistanceConstraint.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <math/MatrixUtils.h>
#include <settings/All.h>

using namespace ausaxs;
using namespace ausaxs::data;
using namespace ausaxs::rigidbody;

TEST_CASE("Backup: Parameters updated in configuration after transformation") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::Linear;

    auto bodies = BodySplitter::split("tests/files/LAR1-2.pdb", {9, 99});
    Rigidbody rigidbody(std::move(bodies));
    rigidbody.molecule.generate_new_hydration();

    unsigned int ibody = 0;
    auto& transformer = rigidbody.transformer;

    auto original_params = rigidbody.conformation->configuration.parameters[ibody];
    Vector3<double> delta_rotation = {1, 1, 1};
    Vector3<double> delta_translation = {1, 1, 1};
    auto delta_params = parameter::RelativeTransformParameters{delta_translation, delta_rotation};
    
    // Expected rotation: R_new = R_delta * R_original (where R_original is identity since original_params.rotation starts at 0)
    auto expected_R = matrix::rotation_matrix(delta_rotation);
    
    // Expected translation: t_new = R_delta * t_original + t_delta
    // Note: original_params.translation is the body's original center of mass, NOT zero!
    auto R_delta = matrix::rotation_matrix(delta_rotation);
    auto expected_translation = R_delta * original_params.translation + delta_translation;

    transformer->apply(std::move(delta_params), ibody);

    // Verify configuration.parameters were updated
    auto& updated_params = rigidbody.conformation->configuration.parameters[ibody];

    // Check that the stored rotation produces the same rotation matrix as the delta
    auto actual_R = matrix::rotation_matrix(updated_params.rotation);
    INFO("Rotation matrices should be equivalent after transformation");
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            REQUIRE_THAT(actual_R(i, j), Catch::Matchers::WithinAbs(expected_R(i, j), 1e-6));
        }
    }

    INFO("Translation should be updated in configuration after transformation");
    REQUIRE_THAT(updated_params.translation.x(), Catch::Matchers::WithinAbs(expected_translation.x(), 1e-6));
    REQUIRE_THAT(updated_params.translation.y(), Catch::Matchers::WithinAbs(expected_translation.y(), 1e-6));
    REQUIRE_THAT(updated_params.translation.z(), Catch::Matchers::WithinAbs(expected_translation.z(), 1e-6));
}

TEST_CASE("Backup: Parameters restored after undo") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::Linear;

    auto bodies = BodySplitter::split("tests/files/LAR1-2.pdb", {9, 99});
    Rigidbody rigidbody(std::move(bodies));
    rigidbody.molecule.generate_new_hydration();
    
    unsigned int ibody = 0;
    auto& transformer = rigidbody.transformer;
    auto& param_gen = rigidbody.parameter_generator;

    // Store original parameters
    auto original_params = rigidbody.conformation->configuration.parameters[ibody];
    auto original_rotation = original_params.rotation;
    auto original_translation = original_params.translation;

    // Apply transformation
    auto new_params = param_gen->next(ibody);
    transformer->apply(std::move(new_params), ibody);
    transformer->undo();

    // Verify parameters were restored
    auto& restored_params = rigidbody.conformation->configuration.parameters[ibody];

    INFO("Rotation should be restored after undo");
    REQUIRE_THAT(restored_params.rotation.x(), Catch::Matchers::WithinAbs(original_rotation.x(), 1e-6));
    REQUIRE_THAT(restored_params.rotation.y(), Catch::Matchers::WithinAbs(original_rotation.y(), 1e-6));
    REQUIRE_THAT(restored_params.rotation.z(), Catch::Matchers::WithinAbs(original_rotation.z(), 1e-6));
    
    INFO("Translation should be restored after undo");
    REQUIRE_THAT(restored_params.translation.x(), Catch::Matchers::WithinAbs(original_translation.x(), 1e-6));
    REQUIRE_THAT(restored_params.translation.y(), Catch::Matchers::WithinAbs(original_translation.y(), 1e-6));
    REQUIRE_THAT(restored_params.translation.z(), Catch::Matchers::WithinAbs(original_translation.z(), 1e-6));
}

TEST_CASE("Backup: Body positions match parameters after transformation") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::Linear;

    auto bodies = BodySplitter::split("tests/files/LAR1-2.pdb", {9, 99});
    Rigidbody rigidbody(std::move(bodies));
    rigidbody.molecule.generate_new_hydration();

    unsigned int ibody = 0;
    auto& transformer = rigidbody.transformer;
    auto& param_gen = rigidbody.parameter_generator;

    // Apply transformation
    auto new_params = param_gen->next(ibody);
    transformer->apply(std::move(new_params), ibody);

    // Reconstruct body from original + parameters
    auto& original_body = rigidbody.conformation->original_conformation[ibody];
    auto& params = rigidbody.conformation->configuration.parameters[ibody];

    // Manually apply transformation to original body
    auto reconstructed_body = original_body;
    reconstructed_body.rotate(matrix::rotation_matrix(params.rotation));
    reconstructed_body.translate(params.translation);

    // Compare with actual body
    auto& actual_body = rigidbody.molecule.get_body(ibody);

    INFO("Body positions should match reconstruction from original + parameters");
    REQUIRE(actual_body.size_atom() == reconstructed_body.size_atom());
    for (unsigned int i = 0; i < actual_body.size_atom(); ++i) {
        auto& actual_atom = actual_body.get_atom(i);
        auto& reconstructed_atom = reconstructed_body.get_atom(i);

        REQUIRE_THAT(actual_atom.coordinates().x(), Catch::Matchers::WithinAbs(reconstructed_atom.coordinates().x(), 1e-4));
        REQUIRE_THAT(actual_atom.coordinates().y(), Catch::Matchers::WithinAbs(reconstructed_atom.coordinates().y(), 1e-4));
        REQUIRE_THAT(actual_atom.coordinates().z(), Catch::Matchers::WithinAbs(reconstructed_atom.coordinates().z(), 1e-4));
    }
}

TEST_CASE("Backup: Constraint-based transforms update all affected body parameters") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::Linear;

    SECTION("SingleTransform updates single body") {
        settings::rigidbody::transform_strategy = settings::rigidbody::TransformationStrategyChoice::SingleTransform;

        auto bodies = BodySplitter::split("tests/files/LAR1-2.pdb", {9, 99});
        Rigidbody rigidbody(std::move(bodies));
        rigidbody.molecule.generate_new_hydration();

        unsigned int ibody = 0;
        auto& transformer = rigidbody.transformer;
        auto& constraint = rigidbody.constraints->distance_constraints_map.at(ibody).at(0).get();

        // Store original parameters
        auto original_params = rigidbody.conformation->configuration.parameters;

        // Apply constraint-based transformation
        Vector3<double> delta_rotation = {1, 1, 1};
        Vector3<double> delta_translation = {1, 1, 1};
        auto delta_params = parameter::RelativeTransformParameters{delta_translation, delta_rotation};
        auto expected_R = matrix::rotation_matrix(delta_rotation);
        
        // For constraint-based transforms, the pivot is the constraining atom position
        // t_new = R_delta * (t_old - pivot) + pivot + t_delta
        auto R_delta = matrix::rotation_matrix(delta_rotation);
        Vector3<double> pivot = constraint.get_atom2().coordinates();
        auto expected_translation = R_delta * (original_params[ibody].translation - pivot) + pivot + delta_translation;
        
        transformer->apply(std::move(delta_params), constraint);

        // Verify the selected body's parameters were updated
        auto& updated_params = rigidbody.conformation->configuration.parameters[ibody];

        INFO("SingleTransform should update only the selected body's parameters");
        for (unsigned int i = 0; i < 3; ++i) {
            if (i == ibody) continue;
            auto& other_params = rigidbody.conformation->configuration.parameters[i];
            REQUIRE((other_params.rotation == original_params[i].rotation && other_params.translation == original_params[i].translation));
        }

        // Check rotation matrix equivalence
        auto actual_R = matrix::rotation_matrix(updated_params.rotation);
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                REQUIRE_THAT(actual_R(i, j), Catch::Matchers::WithinAbs(expected_R(i, j), 1e-6));
            }
        }
        REQUIRE_THAT(updated_params.translation.x(), Catch::Matchers::WithinAbs(expected_translation.x(), 1e-6));
        REQUIRE_THAT(updated_params.translation.y(), Catch::Matchers::WithinAbs(expected_translation.y(), 1e-6));
        REQUIRE_THAT(updated_params.translation.z(), Catch::Matchers::WithinAbs(expected_translation.z(), 1e-6));
    }

    SECTION("RigidTransform updates all connected bodies") {
        settings::rigidbody::transform_strategy = settings::rigidbody::TransformationStrategyChoice::RigidTransform;

        auto bodies = BodySplitter::split("tests/files/LAR1-2.pdb", {9, 99, 199});
        Rigidbody rigidbody(std::move(bodies));
        rigidbody.molecule.generate_new_hydration();

        unsigned int ibody = 1; // Middle body
        auto& transformer = rigidbody.transformer;
        auto& param_gen = rigidbody.parameter_generator;
        auto& constraint = rigidbody.constraints->distance_constraints_map.at(ibody).at(0).get();

        // Store original parameters for all bodies
        std::vector<rigidbody::parameter::BodyTransformParameters> original_params = rigidbody.conformation->configuration.parameters;

        // Apply rigid transformation
        auto new_params = param_gen->next(ibody);
        transformer->apply(std::move(new_params), constraint);

        // At least one body should have updated parameters
        bool any_updated = false;
        for (unsigned int i = 0; i < rigidbody.molecule.size_body(); ++i) {
            auto& updated = rigidbody.conformation->configuration.parameters[i];
            auto& original = original_params[i];

            if (updated.rotation != original.rotation || updated.translation != original.translation) {
                any_updated = true;
                break;
            }
        }

        INFO("RigidTransform should update at least one body's parameters");
        REQUIRE(any_updated);
    }
}

TEST_CASE("Backup: Apply-undo-apply cycle maintains consistency") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::Linear;

    auto bodies = BodySplitter::split("tests/files/LAR1-2.pdb", {9, 99});
    Rigidbody rigidbody(std::move(bodies));
    rigidbody.molecule.generate_new_hydration();

    unsigned int ibody = 0;
    auto& transformer = rigidbody.transformer;
    auto& param_gen = rigidbody.parameter_generator;

    // Store initial state
    auto state0_params = rigidbody.conformation->configuration.parameters[ibody];
    auto state0_body = rigidbody.molecule.get_body(ibody);

    // Apply transformation
    auto new_params = param_gen->next(ibody);
    transformer->apply(std::move(new_params), ibody);
    auto state1_params = rigidbody.conformation->configuration.parameters[ibody];

    // Verify transformation happened
    REQUIRE((state0_params.rotation != state1_params.rotation || state0_params.translation != state1_params.translation));

    // Undo transformation
    transformer->undo();
    auto restored_params = rigidbody.conformation->configuration.parameters[ibody];  // Make a copy, not a reference
    auto& restored_body = rigidbody.molecule.get_body(ibody);

    INFO("Undo should restore to previous state");
    REQUIRE_THAT(restored_params.rotation.x(), Catch::Matchers::WithinAbs(state0_params.rotation.x(), 1e-6));
    REQUIRE_THAT(restored_params.rotation.y(), Catch::Matchers::WithinAbs(state0_params.rotation.y(), 1e-6));
    REQUIRE_THAT(restored_params.rotation.z(), Catch::Matchers::WithinAbs(state0_params.rotation.z(), 1e-6));
    REQUIRE_THAT(restored_params.translation.x(), Catch::Matchers::WithinAbs(state0_params.translation.x(), 1e-6));
    REQUIRE_THAT(restored_params.translation.y(), Catch::Matchers::WithinAbs(state0_params.translation.y(), 1e-6));
    REQUIRE_THAT(restored_params.translation.z(), Catch::Matchers::WithinAbs(state0_params.translation.z(), 1e-6));

    // Apply a new transformation
    auto new_params2 = param_gen->next(ibody);
    transformer->apply(std::move(new_params2), ibody);
    auto state2_params = rigidbody.conformation->configuration.parameters[ibody];

    INFO("After undo and new apply, we should be in a different state");
    REQUIRE((restored_params.rotation != state2_params.rotation || restored_params.translation != state2_params.translation));

    // Undo again - should go back to state after first undo
    transformer->undo();
    auto& final_params = rigidbody.conformation->configuration.parameters[ibody];

    INFO("Second undo should restore to state before second apply (which was state0)");
    REQUIRE_THAT(final_params.rotation.x(), Catch::Matchers::WithinAbs(state0_params.rotation.x(), 1e-6));
    REQUIRE_THAT(final_params.rotation.y(), Catch::Matchers::WithinAbs(state0_params.rotation.y(), 1e-6));
    REQUIRE_THAT(final_params.rotation.z(), Catch::Matchers::WithinAbs(state0_params.rotation.z(), 1e-6));
    REQUIRE_THAT(final_params.translation.x(), Catch::Matchers::WithinAbs(state0_params.translation.x(), 1e-6));
    REQUIRE_THAT(final_params.translation.y(), Catch::Matchers::WithinAbs(state0_params.translation.y(), 1e-6));
    REQUIRE_THAT(final_params.translation.z(), Catch::Matchers::WithinAbs(state0_params.translation.z(), 1e-6));
}
