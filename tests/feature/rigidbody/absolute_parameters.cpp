#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <rigidbody/Rigidbody.h>
#include <rigidbody/BodySplitter.h>
#include <rigidbody/controller/SimpleController.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/detail/Conformation.h>
#include <rigidbody/detail/Configuration.h>
#include <rigidbody/parameters/BodyTransformParameters.h>
#include <rigidbody/parameters/ParameterGenerationStrategy.h>
#include <rigidbody/transform/TransformStrategy.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <math/MatrixUtils.h>
#include <settings/All.h>
#include <io/ExistingFile.h>

using namespace ausaxs;
using namespace ausaxs::data;
using namespace ausaxs::rigidbody;

/**
 * @brief Verify that applying configuration.parameters to original_conformation reproduces the current body state.
 */
void verify_configuration_consistency(const Rigidbody& rigidbody) {
    for (size_t ibody = 0; ibody < rigidbody.molecule.size_body(); ++ibody) {
        const auto& current_body = rigidbody.molecule.get_body(ibody);
        const auto& original_body = rigidbody.conformation->original_conformation[ibody];
        const auto& params = rigidbody.conformation->configuration.parameters[ibody];

        // original_conformation should be centered at origin
        auto original_cm = original_body.get_cm(false);
        REQUIRE_THAT(original_cm.x(), Catch::Matchers::WithinAbs(0.0, 1e-6));
        REQUIRE_THAT(original_cm.y(), Catch::Matchers::WithinAbs(0.0, 1e-6));
        REQUIRE_THAT(original_cm.z(), Catch::Matchers::WithinAbs(0.0, 1e-6));

        // reconstruct the body from original_conformation + parameters
        Body reconstructed = original_body; // copy
        reconstructed.rotate(matrix::rotation_matrix(params.rotation));
        reconstructed.translate(params.translation);

        // reconstructed body should match current body
        auto current_cm = current_body.get_cm(false);
        auto reconstructed_cm = reconstructed.get_cm(false);

        INFO("Body " << ibody << " center of mass mismatch");
        REQUIRE_THAT(reconstructed_cm.x(), Catch::Matchers::WithinAbs(current_cm.x(), 1e-3));
        REQUIRE_THAT(reconstructed_cm.y(), Catch::Matchers::WithinAbs(current_cm.y(), 1e-3));
        REQUIRE_THAT(reconstructed_cm.z(), Catch::Matchers::WithinAbs(current_cm.z(), 1e-3));

        // individual atoms should match
        REQUIRE(reconstructed.size_atom() == current_body.size_atom());
        for (size_t iatom = 0; iatom < current_body.size_atom(); ++iatom) {
            auto current_pos = current_body.get_atom(iatom).coordinates();
            auto reconstructed_pos = reconstructed.get_atom(iatom).coordinates();

            INFO("Body " << ibody << ", atom " << iatom << " position mismatch");
            REQUIRE_THAT(reconstructed_pos.x(), Catch::Matchers::WithinAbs(current_pos.x(), 1e-3));
            REQUIRE_THAT(reconstructed_pos.y(), Catch::Matchers::WithinAbs(current_pos.y(), 1e-3));
            REQUIRE_THAT(reconstructed_pos.z(), Catch::Matchers::WithinAbs(current_pos.z(), 1e-3));
        }
    }
}

TEST_CASE("AbsoluteParameters: Initial configuration consistency") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::Linear;

    SECTION("simple bodies") {
        // create a simple multi-body system
        AtomFF a1({0, 0, 0}, form_factor::form_factor_t::C);
        AtomFF a2({1, 0, 0}, form_factor::form_factor_t::C);
        AtomFF a3({0, 1, 0}, form_factor::form_factor_t::C);
        AtomFF a4({0, 0, 1}, form_factor::form_factor_t::C);

        Body b1 = Body(std::vector<AtomFF>{a1});
        Body b2 = Body(std::vector<AtomFF>{a2});
        Body b3 = Body(std::vector<AtomFF>{a3});
        Body b4 = Body(std::vector<AtomFF>{a4});

        Rigidbody rigidbody(Molecule{std::vector<Body>{b1, b2, b3, b4}});

        // initial configuration should satisfy the invariant
        verify_configuration_consistency(rigidbody);
    }

    SECTION("real protein structure") {
        Rigidbody rigidbody = BodySplitter::split("tests/files/LAR1-2.pdb", {9, 99});

        // initial configuration should satisfy the invariant
        verify_configuration_consistency(rigidbody);
    }
}

TEST_CASE("AbsoluteParameters: Free body transformations preserve consistency") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::None;
    settings::rigidbody::iterations = 10;
    settings::grid::min_bins = 100;

    // simple system with no constraints
    AtomFF a1({0, 0, 0}, form_factor::form_factor_t::C);
    AtomFF a2({1, 0, 0}, form_factor::form_factor_t::C);
    AtomFF a3({0, 1, 0}, form_factor::form_factor_t::C);

    Body b1 = Body(std::vector<AtomFF>{a1});
    Body b2 = Body(std::vector<AtomFF>{a2});
    Body b3 = Body(std::vector<AtomFF>{a3});

    Rigidbody rigidbody(Molecule{std::vector<Body>{b1, b2, b3}});

    // apply some transformations using the base class apply(param, ibody)
    auto& transformer = rigidbody.transformer;
    auto& param_gen = rigidbody.parameter_generator;
    for (int iter = 0; iter < 5; ++iter) {
        for (size_t ibody = 0; ibody < rigidbody.molecule.size_body(); ++ibody) {
            auto params = param_gen->next(ibody);
            transformer->apply(std::move(params), ibody);

            INFO("After iteration " << iter << ", body " << ibody);
            verify_configuration_consistency(rigidbody);
        }
    }
}

TEST_CASE("AbsoluteParameters: Constraint-based transformations preserve consistency") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::Linear;
    settings::rigidbody::transform_strategy = settings::rigidbody::TransformationStrategyChoice::RigidTransform;
    settings::rigidbody::iterations = 10;

    SECTION("SingleTransform strategy") {
        settings::rigidbody::transform_strategy = settings::rigidbody::TransformationStrategyChoice::SingleTransform;

        Rigidbody rigidbody = BodySplitter::split("tests/files/LAR1-2.pdb", {9, 99});
        verify_configuration_consistency(rigidbody);

        // apply transformations via constraints
        auto& transformer = rigidbody.transformer;
        auto& param_gen = rigidbody.parameter_generator;

        for (int iter = 0; iter < 5; ++iter) {
            // transform the first constrained body
            unsigned int ibody = 0;
            auto& constraint = rigidbody.constraints->distance_constraints_map.at(ibody).at(0).get();
            auto params = param_gen->next(ibody);
            
            transformer->apply(std::move(params), constraint);
            
            INFO("After iteration " << iter << " with SingleTransform");
            verify_configuration_consistency(rigidbody);
        }
    }

    SECTION("RigidTransform strategy") {
        settings::rigidbody::transform_strategy = settings::rigidbody::TransformationStrategyChoice::RigidTransform;
        
        Rigidbody rigidbody = BodySplitter::split("tests/files/LAR1-2.pdb", {9, 99});
        verify_configuration_consistency(rigidbody);

        // apply transformations via constraints
        auto& transformer = rigidbody.transformer;
        auto& param_gen = rigidbody.parameter_generator;

        for (int iter = 0; iter < 5; ++iter) {
            // transform the first constrained body
            unsigned int ibody = 0;
            auto& constraint = rigidbody.constraints->distance_constraints_map.at(ibody).at(0).get();
            auto params = param_gen->next(ibody);
            
            transformer->apply(std::move(params), constraint);
            
            INFO("After iteration " << iter << " with RigidTransform");
            verify_configuration_consistency(rigidbody);
        }
    }
}

TEST_CASE("AbsoluteParameters: Full optimization run preserves consistency") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::Linear;
    settings::rigidbody::iterations = 20;

    SECTION("with SingleTransform") {
        settings::rigidbody::transform_strategy = settings::rigidbody::TransformationStrategyChoice::SingleTransform;
        
        Rigidbody rigidbody = BodySplitter::split("tests/files/LAR1-2.pdb", {9, 99});
        auto controller = rigidbody.controller.get();
        controller->setup(io::ExistingFile("tests/files/LAR1-2.dat"));

        // run several optimization steps
        for (int i = 0; i < 10; ++i) {
            controller->run_step();
            
            INFO("After optimization step " << i << " with SingleTransform");
            verify_configuration_consistency(rigidbody);
        }
    }

    SECTION("with RigidTransform") {
        settings::rigidbody::transform_strategy = settings::rigidbody::TransformationStrategyChoice::RigidTransform;
        
        Rigidbody rigidbody = BodySplitter::split("tests/files/LAR1-2.pdb", {9, 99});
        auto controller = rigidbody.controller.get();
        controller->setup(io::ExistingFile("tests/files/LAR1-2.dat"));

        // run several optimization steps
        for (int i = 0; i < 10; ++i) {
            controller->run_step();
            
            INFO("After optimization step " << i << " with RigidTransform");
            verify_configuration_consistency(rigidbody);
        }
    }
}

TEST_CASE("AbsoluteParameters: Undo restores configuration.parameters") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::None;

    AtomFF a1({0, 0, 0}, form_factor::form_factor_t::C);
    AtomFF a2({1, 0, 0}, form_factor::form_factor_t::C);
    Body b1 = Body(std::vector<AtomFF>{a1});
    Body b2 = Body(std::vector<AtomFF>{a2});

    Rigidbody rigidbody(Molecule{std::vector<Body>{b1, b2}});

    // store original parameters
    auto original_params = rigidbody.conformation->configuration.parameters[0];

    // apply a transformation
    auto& transformer = rigidbody.transformer;
    auto& param_gen = rigidbody.parameter_generator;
    auto new_params = param_gen->next(0);

    // store the new parameters for comparison
    auto expected_new_translation = new_params.translation;
    auto expected_new_rotation = new_params.rotation;

    transformer->apply(std::move(new_params), 0);
    
    // verify parameters were updated
    auto& updated_params = rigidbody.conformation->configuration.parameters[0];
    REQUIRE(updated_params.translation.x() == expected_new_translation.x());
    REQUIRE(updated_params.translation.y() == expected_new_translation.y());
    REQUIRE(updated_params.translation.z() == expected_new_translation.z());

    // undo the transformation
    transformer->undo();

    // verify parameters were restored
    auto& restored_params = rigidbody.conformation->configuration.parameters[0];
    REQUIRE_THAT(restored_params.translation.x(), Catch::Matchers::WithinAbs(original_params.translation.x(), 1e-6));
    REQUIRE_THAT(restored_params.translation.y(), Catch::Matchers::WithinAbs(original_params.translation.y(), 1e-6));
    REQUIRE_THAT(restored_params.translation.z(), Catch::Matchers::WithinAbs(original_params.translation.z(), 1e-6));
    REQUIRE_THAT(restored_params.rotation.x(), Catch::Matchers::WithinAbs(original_params.rotation.x(), 1e-6));
    REQUIRE_THAT(restored_params.rotation.y(), Catch::Matchers::WithinAbs(original_params.rotation.y(), 1e-6));
    REQUIRE_THAT(restored_params.rotation.z(), Catch::Matchers::WithinAbs(original_params.rotation.z(), 1e-6));

    // and verify consistency is maintained
    verify_configuration_consistency(rigidbody);
}
