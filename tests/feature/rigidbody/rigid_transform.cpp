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
    rigidbody.molecule.generate_new_hydration();

    auto& transformer = rigidbody.transformer;
    auto& param_gen = rigidbody.parameter_generator;
    
    // Transform via a constraint - RigidTransform will choose which side to transform
    unsigned int ibody = 1;
    auto& constraint = rigidbody.constraints->distance_constraints_map.at(ibody).at(0).get();

    // Store initial parameters for all bodies
    std::vector<rigidbody::parameter::BodyTransformParameters> initial_params = rigidbody.conformation->configuration.parameters;

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

        INFO("Body " << i << " should be reconstructible from parameters");
        REQUIRE_THAT(reconstructed_cm.x(), Catch::Matchers::WithinAbs(current_cm.x(), 1e-3));
        REQUIRE_THAT(reconstructed_cm.y(), Catch::Matchers::WithinAbs(current_cm.y(), 1e-3));
        REQUIRE_THAT(reconstructed_cm.z(), Catch::Matchers::WithinAbs(current_cm.z(), 1e-3));
    }
}

TEST_CASE("RigidTransform: Distance constraints preserved") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::Linear;
    settings::rigidbody::transform_strategy = settings::rigidbody::TransformationStrategyChoice::RigidTransform;

    auto bodies = BodySplitter::split("tests/files/LAR1-2.pdb", {9, 99, 199});
    Rigidbody rigidbody(std::move(bodies));
    rigidbody.molecule.generate_new_hydration();

    auto& transformer = rigidbody.transformer;
    auto& param_gen = rigidbody.parameter_generator;

    // Record initial constraint distances
    std::vector<double> initial_distances;
    for (const auto& constraint : rigidbody.constraints->distance_constraints) {
        auto dist = (constraint.get_atom1().coordinates() - constraint.get_atom2().coordinates()).norm();
        initial_distances.push_back(dist);
    }

    // Apply several transformations
    for (int iter = 0; iter < 5; ++iter) {
        for (size_t ibody = 0; ibody < rigidbody.molecule.size_body(); ++ibody) {
            if (rigidbody.constraints->distance_constraints_map.at(ibody).empty()) continue;
            
            auto& constraint = rigidbody.constraints->distance_constraints_map.at(ibody).at(0).get();
            auto params = param_gen->next(ibody);
            transformer->apply(std::move(params), constraint);

            // Verify all constraint distances remain approximately the same
            for (size_t i = 0; i < rigidbody.constraints->distance_constraints.size(); ++i) {
                auto& c = rigidbody.constraints->distance_constraints[i];
                double current_distance = (c.get_atom1().coordinates() - c.get_atom2().coordinates()).norm();
                INFO("After iteration " << iter << ", body " << ibody << ", constraint " << i);
                REQUIRE_THAT(current_distance, Catch::Matchers::WithinAbs(initial_distances[i], 0.1));
            }
        }
    }
}

TEST_CASE("RigidTransform: Orbital motion correctness") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::None;
    settings::molecule::center = false;
    settings::grid::scaling = 2;

    // Create simple test bodies with known geometry - bodies at x = -5, 0, +5
    AtomFF a1({-5, 0, 0}, form_factor::form_factor_t::C);
    AtomFF a2({0, 0, 0}, form_factor::form_factor_t::C);
    AtomFF a3({5, 0, 0}, form_factor::form_factor_t::C);
    
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
    rigidbody.molecule.generate_new_hydration();

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
