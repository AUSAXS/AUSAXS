#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <rigidbody/Rigidbody.h>
#include <rigidbody/BodySplitter.h>
#include <rigidbody/transform/RigidTransform.h>
#include <rigidbody/constraints/DistanceConstraint.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/parameters/BodyTransformParameters.h>
#include <rigidbody/parameters/ParameterGenerationStrategy.h>
#include <rigidbody/detail/Conformation.h>
#include <data/Body.h>
#include <data/Molecule.h>
#include <math/MatrixUtils.h>
#include <settings/All.h>

#include <numbers>

using namespace ausaxs;
using namespace ausaxs::data;
using namespace ausaxs::rigidbody;

TEST_CASE("RigidTransform: Secondary body parameter updates") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::Linear;
    settings::rigidbody::transform_strategy = settings::rigidbody::TransformationStrategyChoice::RigidTransform;

    auto bodies = BodySplitter::split("tests/files/LAR1-2.pdb", {9, 99, 199});
    Rigidbody rigidbody(std::move(bodies));
    // Note: Do not call generate_new_hydration() here - it changes CM after Conformation is created

    auto& transformer = rigidbody.transformer;
    auto& param_gen = rigidbody.parameter_generator;
    
    // Transform via a constraint - RigidTransform will choose which side to transform
    unsigned int ibody = 1;
    auto& constraint = rigidbody.constraints->distance_constraints_map.at(ibody).at(0).get();

    // Store initial parameters for all bodies
    std::vector<rigidbody::parameter::BodyTransformParameters> initial_params = rigidbody.conformation->configuration.parameters;
    
    // Also store initial body CMs
    std::vector<Vector3<double>> initial_cms;
    for (unsigned int i = 0; i < rigidbody.molecule.size_body(); ++i) {
        initial_cms.push_back(rigidbody.molecule.get_body(i).get_cm());
    }

    // Apply transformation
    auto new_params = param_gen->next(ibody);
    transformer->apply(std::move(new_params), constraint);

    // At least one body should have updated parameters (the rigid group that was transformed)
    int bodies_updated = 0;
    for (size_t i = 0; i < rigidbody.molecule.size_body(); ++i) {
        auto& updated = rigidbody.conformation->configuration.parameters[i];
        auto& original = initial_params[i];

        if (updated.rotation != original.rotation || updated.translation != original.translation) {
            bodies_updated++;
        }
    }

    INFO("RigidTransform should update at least one body's parameters");
    REQUIRE(bodies_updated >= 1);

    // Verify parameters can reconstruct the current state for ALL bodies
    for (size_t i = 0; i < rigidbody.molecule.size_body(); ++i) {
        auto& current_body = rigidbody.molecule.get_body(i);
        auto& params = rigidbody.conformation->configuration.parameters[i];
        auto& original = rigidbody.conformation->original_conformation[i];

        // Reconstruct from original + parameters
        Body reconstructed = original;
        reconstructed.rotate(matrix::rotation_matrix(params.rotation));
        reconstructed.translate(params.translation);

        auto current_cm = current_body.get_cm();
        auto reconstructed_cm = reconstructed.get_cm();

        INFO("Body " << i << " reconstruction should match current position");
        REQUIRE_THAT(reconstructed_cm.x(), Catch::Matchers::WithinAbs(current_cm.x(), 1e-3));
        REQUIRE_THAT(reconstructed_cm.y(), Catch::Matchers::WithinAbs(current_cm.y(), 1e-3));
        REQUIRE_THAT(reconstructed_cm.z(), Catch::Matchers::WithinAbs(current_cm.z(), 1e-3));
    }
}

TEST_CASE("RigidTransform: Internal constraints within group preserved") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::Linear;
    settings::rigidbody::transform_strategy = settings::rigidbody::TransformationStrategyChoice::RigidTransform;

    auto bodies = BodySplitter::split("tests/files/LAR1-2.pdb", {9, 99, 199});
    Rigidbody rigidbody(std::move(bodies));
    rigidbody.molecule.generate_new_hydration();

    auto& transformer = rigidbody.transformer;
    auto& param_gen = rigidbody.parameter_generator;

    // For a chain A-B-C with constraints A-B and B-C:
    // When transforming at constraint A-B:
    //   - If body A is moved, constraint A-B distance can change (it's the hinge)
    //   - Constraint B-C should be preserved (B and C stay together)
    // When transforming at constraint B-C:
    //   - Constraint A-B should be preserved
    //   - Constraint B-C distance can change

    // Record initial constraint distances
    std::vector<double> initial_distances;
    for (const auto& constraint : rigidbody.constraints->distance_constraints) {
        auto dist = (constraint.get_atom1().coordinates() - constraint.get_atom2().coordinates()).norm();
        initial_distances.push_back(dist);
    }

    // Transform using constraint 0 (between body 0 and 1)
    // This should preserve constraint 1 (between body 1 and 2)
    auto& constraint0 = rigidbody.constraints->distance_constraints[0];
    auto params = param_gen->next(0);
    transformer->apply(std::move(params), constraint0);

    // Check constraint 1 (not the hinge) is preserved
    auto& c1 = rigidbody.constraints->distance_constraints[1];
    double new_distance_1 = (c1.get_atom1().coordinates() - c1.get_atom2().coordinates()).norm();
    INFO("Constraint 1 (internal to non-moving group) should be preserved");
    REQUIRE_THAT(new_distance_1, Catch::Matchers::WithinAbs(initial_distances[1], 0.1));
}

TEST_CASE("RigidTransform: Orbital motion correctness") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::None;
    settings::molecule::center = false;
    settings::grid::scaling = 2;

    // Create simple test bodies with known geometry - bodies at x = -3, 0, +3
    // Atoms must be within 4 angstroms for valid distance constraints
    AtomFF a1({-3, 0, 0}, form_factor::form_factor_t::C);
    AtomFF a2({0, 0, 0}, form_factor::form_factor_t::C);
    AtomFF a3({3, 0, 0}, form_factor::form_factor_t::C);
    
    Body b1 = Body(std::vector<AtomFF>{a1});
    Body b2 = Body(std::vector<AtomFF>{a2});
    Body b3 = Body(std::vector<AtomFF>{a3});

    Rigidbody rigidbody(Molecule{std::vector<Body>{b1, b2, b3}});
    
    // Manually create constraints: 0 - 1 - 2 (only 0-1 constraint so transforming it affects only body 0)
    rigidbody.constraints->add_constraint(
        rigidbody::constraints::DistanceConstraint(&rigidbody.molecule, 0, 1, 0, 0)
    );

    auto& transformer = rigidbody.transformer;

    // Record initial distances
    auto initial_cm_0 = rigidbody.molecule.get_body(0).get_cm();
    auto initial_cm_1 = rigidbody.molecule.get_body(1).get_cm();
    double initial_dist = (initial_cm_0 - initial_cm_1).norm();

    // Rotate body 0 by 90 degrees around Z axis (body 1 is at origin)
    auto& constraint = rigidbody.constraints->distance_constraints[0];
    transformer->apply({{0, 0, 0}, {0, 0, std::numbers::pi/2}}, constraint);

    // Verify distance is preserved (rigid relationship maintained)
    auto final_cm_0 = rigidbody.molecule.get_body(0).get_cm();
    auto final_cm_1 = rigidbody.molecule.get_body(1).get_cm();
    double final_dist = (final_cm_0 - final_cm_1).norm();

    INFO("Rigid relationship should be preserved during rotation");
    REQUIRE_THAT(final_dist, Catch::Matchers::WithinAbs(initial_dist, 1e-6));

    // Verify transformed body had its parameters updated and can be reconstructed
    auto& params_0 = rigidbody.conformation->configuration.parameters[0];
    auto& original_0 = rigidbody.conformation->original_conformation[0];
    
    Body reconstructed = original_0;
    reconstructed.rotate(matrix::rotation_matrix(params_0.rotation));
    reconstructed.translate(params_0.translation);
    
    auto current_cm = rigidbody.molecule.get_body(0).get_cm();
    auto reconstructed_cm = reconstructed.get_cm();
    
    INFO("Transformed body orbital motion should be captured in parameters");
    REQUIRE_THAT(reconstructed_cm.x(), Catch::Matchers::WithinAbs(current_cm.x(), 1e-3));
    REQUIRE_THAT(reconstructed_cm.y(), Catch::Matchers::WithinAbs(current_cm.y(), 1e-3));
    REQUIRE_THAT(reconstructed_cm.z(), Catch::Matchers::WithinAbs(current_cm.z(), 1e-3));
}

TEST_CASE("RigidTransform: Multi-step transformation consistency") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::Linear;
    settings::rigidbody::transform_strategy = settings::rigidbody::TransformationStrategyChoice::RigidTransform;

    auto bodies = BodySplitter::split("tests/files/LAR1-2.pdb", {9, 99});
    Rigidbody rigidbody(std::move(bodies));
    // Note: Do not call generate_new_hydration() here - it changes CM after Conformation is created

    auto& transformer = rigidbody.transformer;
    auto& param_gen = rigidbody.parameter_generator;

    // Apply multiple transformations and verify consistency is maintained throughout
    for (int iter = 0; iter < 10; ++iter) {
        unsigned int ibody = iter % rigidbody.molecule.size_body();
        if (rigidbody.constraints->distance_constraints_map.at(ibody).empty()) continue;

        auto& constraint = rigidbody.constraints->distance_constraints_map.at(ibody).at(0).get();
        auto params = param_gen->next(ibody);
        transformer->apply(std::move(params), constraint);

        // Verify all bodies can be reconstructed from parameters
        for (size_t i = 0; i < rigidbody.molecule.size_body(); ++i) {
            auto& current_body = rigidbody.molecule.get_body(i);
            auto& body_params = rigidbody.conformation->configuration.parameters[i];
            auto& original = rigidbody.conformation->original_conformation[i];

            Body reconstructed = original;
            reconstructed.rotate(matrix::rotation_matrix(body_params.rotation));
            reconstructed.translate(body_params.translation);

            auto current_cm = current_body.get_cm();
            auto reconstructed_cm = reconstructed.get_cm();

            INFO("After iteration " << iter << ", body " << i);
            REQUIRE_THAT(reconstructed_cm.x(), Catch::Matchers::WithinAbs(current_cm.x(), 1e-3));
            REQUIRE_THAT(reconstructed_cm.y(), Catch::Matchers::WithinAbs(current_cm.y(), 1e-3));
            REQUIRE_THAT(reconstructed_cm.z(), Catch::Matchers::WithinAbs(current_cm.z(), 1e-3));
        }
    }
}
