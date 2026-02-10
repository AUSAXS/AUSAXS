#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <rigidbody/Rigidbody.h>
#include <rigidbody/BodySplitter.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/detail/SystemSpecification.h>
#include <rigidbody/detail/MoleculeTransformParametersAbsolute.h>
#include <rigidbody/parameters/BodyTransformParametersAbsolute.h>
#include <rigidbody/parameters/OptimizableSymmetryStorage.h>
#include <rigidbody/parameters/ParameterGenerationStrategy.h>
#include <rigidbody/transform/TransformStrategy.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <data/symmetry/BodySymmetryFacade.h>
#include <settings/All.h>
#include <io/ExistingFile.h>

using namespace ausaxs;
using namespace ausaxs::data;
using namespace ausaxs::rigidbody;

TEST_CASE("SymmetryBackup: Symmetry structure preserved in original_conformation") {
    settings::general::verbose = false;
    settings::grid::min_bins = 250;
    settings::molecule::implicit_hydrogens = false;

    SECTION("via direct construction") {
        auto bodies = BodySplitter::split("tests/files/LAR1-2.pdb", {9, 99});
        bodies.get_body(0).symmetry().add(symmetry::type::p2);
        
        Rigidbody rigidbody(std::move(bodies));
        rigidbody.molecule.generate_new_hydration();
        REQUIRE(rigidbody.molecule.get_body(0).size_symmetry() == 1);

        // Verify original_conformation also has matching symmetry structure
        auto& original_body = rigidbody.conformation->initial_conformation[0];
        INFO("original_conformation must have same symmetry structure as molecule body");
        REQUIRE(original_body.size_symmetry() == 1);

        // Verify both use OptimizableSymmetryStorage
        auto* mol_storage = dynamic_cast<symmetry::OptimizableSymmetryStorage*>(
            rigidbody.molecule.get_body(0).symmetry().get_obj()
        );
        auto* orig_storage = dynamic_cast<symmetry::OptimizableSymmetryStorage*>(
            original_body.symmetry().get_obj()
        );

        REQUIRE(mol_storage != nullptr);
        REQUIRE(orig_storage != nullptr);

        // Verify configuration.parameters has symmetry_pars initialized
        REQUIRE(rigidbody.conformation->absolute_parameters.parameters[0].symmetry_pars.size() == 1);
    }
}

TEST_CASE("SymmetryBackup: Symmetry parameters backed up and restored on undo") {
    settings::general::verbose = false;
    settings::grid::min_bins = 250;
    settings::molecule::implicit_hydrogens = false;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::Linear;

    // Create a rigidbody with symmetry
    auto bodies = BodySplitter::split("tests/files/LAR1-2.pdb", {9, 99});
    bodies.get_body(0).symmetry().add(symmetry::type::p2);
    
    Rigidbody rigidbody(std::move(bodies));
    rigidbody.molecule.generate_new_hydration();
    unsigned int ibody = 0;

    // Store original symmetry parameters
    auto original_sym_pars = rigidbody.conformation->absolute_parameters.parameters[ibody].symmetry_pars;
    REQUIRE(original_sym_pars.size() == 1);

    auto original_sym_translation = original_sym_pars[0].initial_relation.translation;
    auto original_sym_orientation = original_sym_pars[0].initial_relation.orientation;

    // Apply a transformation that modifies symmetry parameters
    auto& transformer = rigidbody.transformer;
    auto& param_gen = rigidbody.parameter_generator;

    // Generate new parameters (which should include modified symmetry parameters)
    auto new_params = param_gen->next(ibody);

    // Store the new symmetry parameters for verification
    auto expected_new_sym_pars = new_params.symmetry_pars;
    REQUIRE(expected_new_sym_pars.has_value());
    REQUIRE(expected_new_sym_pars.value().size() == 1);

    // Apply the transformation
    transformer->apply(std::move(new_params), ibody);

    // Verify symmetry parameters were updated in configuration
    auto& updated_sym_pars = rigidbody.conformation->absolute_parameters.parameters[ibody].symmetry_pars;
    REQUIRE(updated_sym_pars.size() == 1);

    INFO("Symmetry parameters should be updated after transformation");
    REQUIRE(updated_sym_pars[0].initial_relation.translation.x() == expected_new_sym_pars.value()[0].initial_relation.translation.x());
    REQUIRE(updated_sym_pars[0].initial_relation.translation.y() == expected_new_sym_pars.value()[0].initial_relation.translation.y());
    REQUIRE(updated_sym_pars[0].initial_relation.translation.z() == expected_new_sym_pars.value()[0].initial_relation.translation.z());

    // Undo the transformation
    transformer->undo();

    // Verify symmetry parameters were restored
    auto& restored_sym_pars = rigidbody.conformation->absolute_parameters.parameters[ibody].symmetry_pars;
    REQUIRE(restored_sym_pars.size() == 1);

    INFO("Symmetry parameters should be restored after undo");
    REQUIRE_THAT(restored_sym_pars[0].initial_relation.translation.x(), Catch::Matchers::WithinAbs(original_sym_translation.x(), 1e-6));
    REQUIRE_THAT(restored_sym_pars[0].initial_relation.translation.y(), Catch::Matchers::WithinAbs(original_sym_translation.y(), 1e-6));
    REQUIRE_THAT(restored_sym_pars[0].initial_relation.translation.z(), Catch::Matchers::WithinAbs(original_sym_translation.z(), 1e-6));
    REQUIRE_THAT(restored_sym_pars[0].initial_relation.orientation.x(), Catch::Matchers::WithinAbs(original_sym_orientation.x(), 1e-6));
    REQUIRE_THAT(restored_sym_pars[0].initial_relation.orientation.y(), Catch::Matchers::WithinAbs(original_sym_orientation.y(), 1e-6));
    REQUIRE_THAT(restored_sym_pars[0].initial_relation.orientation.z(), Catch::Matchers::WithinAbs(original_sym_orientation.z(), 1e-6));
}

TEST_CASE("SymmetryBackup: Body symmetry storage preserved through transformations") {
    settings::general::verbose = false;
    settings::grid::min_bins = 250;
    settings::molecule::implicit_hydrogens = false;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::Linear;

    auto bodies = BodySplitter::split("tests/files/LAR1-2.pdb", {9, 99});
    bodies.get_body(0).symmetry().add(symmetry::type::p2);
    
    Rigidbody rigidbody(std::move(bodies));
    rigidbody.molecule.generate_new_hydration();
    unsigned int ibody = 0;

    // Verify initial state has OptimizableSymmetryStorage
    auto* initial_storage = dynamic_cast<symmetry::OptimizableSymmetryStorage*>(
        rigidbody.molecule.get_body(ibody).symmetry().get_obj()
    );
    REQUIRE(initial_storage != nullptr);
    
    // Verify configuration was properly initialized with symmetry parameters
    REQUIRE(rigidbody.conformation->absolute_parameters.parameters[ibody].symmetry_pars.size() == 1);

    // Apply a transformation
    auto& transformer = rigidbody.transformer;
    auto& param_gen = rigidbody.parameter_generator;
    auto new_params = param_gen->next(ibody);
    transformer->apply(std::move(new_params), ibody);

    // Verify the body still has OptimizableSymmetryStorage after transformation
    auto* after_storage = dynamic_cast<symmetry::OptimizableSymmetryStorage*>(
        rigidbody.molecule.get_body(ibody).symmetry().get_obj()
    );
    INFO("Body should retain OptimizableSymmetryStorage after transformation");
    REQUIRE(after_storage != nullptr);
    REQUIRE(rigidbody.molecule.get_body(ibody).size_symmetry() == 1);

    // Undo the transformation
    transformer->undo();

    // Verify the body still has OptimizableSymmetryStorage after undo
    auto* undo_storage = dynamic_cast<symmetry::OptimizableSymmetryStorage*>(
        rigidbody.molecule.get_body(ibody).symmetry().get_obj()
    );
    INFO("Body should retain OptimizableSymmetryStorage after undo");
    REQUIRE(undo_storage != nullptr);
    REQUIRE(rigidbody.molecule.get_body(ibody).size_symmetry() == 1);
}

TEST_CASE("SymmetryBackup: Constraint-based transforms preserve symmetries") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::Linear;

    SECTION("SingleTransform") {
        settings::grid::min_bins = 250;
        settings::rigidbody::transform_strategy = settings::rigidbody::TransformationStrategyChoice::SingleTransform;

        auto bodies = BodySplitter::split("tests/files/LAR1-2.pdb", {9, 99});
        bodies.get_body(0).symmetry().add(symmetry::type::p2);
        
        Rigidbody rigidbody(std::move(bodies));
        rigidbody.molecule.generate_new_hydration();
        unsigned int ibody = 0;

        // Store original symmetry state
        auto original_size = rigidbody.molecule.get_body(ibody).size_symmetry();
        auto original_sym_pars = rigidbody.conformation->absolute_parameters.parameters[ibody].symmetry_pars;

        // Apply constraint-based transformation
        auto& transformer = rigidbody.transformer;
        auto& param_gen = rigidbody.parameter_generator;
        auto constraint = rigidbody.constraints->get_body_constraints(ibody).at(0);

        auto new_params = param_gen->next(ibody);
        transformer->apply(std::move(new_params), constraint);

        // Verify symmetry count is preserved
        INFO("Symmetry count should be preserved after constraint transformation");
        REQUIRE(rigidbody.molecule.get_body(ibody).size_symmetry() == original_size);

        // Verify OptimizableSymmetryStorage is still there
        auto* storage = dynamic_cast<symmetry::OptimizableSymmetryStorage*>(
            rigidbody.molecule.get_body(ibody).symmetry().get_obj()
        );
        REQUIRE(storage != nullptr);

        // Undo and verify restoration
        transformer->undo();

        REQUIRE(rigidbody.molecule.get_body(ibody).size_symmetry() == original_size);
        auto& restored_sym_pars = rigidbody.conformation->absolute_parameters.parameters[ibody].symmetry_pars;
        REQUIRE(restored_sym_pars.size() == original_sym_pars.size());
    }

    SECTION("RigidTransform") {
        settings::grid::min_bins = 250;
        settings::rigidbody::transform_strategy = settings::rigidbody::TransformationStrategyChoice::RigidTransform;

        auto bodies = BodySplitter::split("tests/files/LAR1-2.pdb", {9, 99});
        bodies.get_body(0).symmetry().add(symmetry::type::p2);
        
        Rigidbody rigidbody(std::move(bodies));
        rigidbody.molecule.generate_new_hydration();
        unsigned int ibody = 0;

        // Store original symmetry state
        auto original_size = rigidbody.molecule.get_body(ibody).size_symmetry();
        auto original_sym_pars = rigidbody.conformation->absolute_parameters.parameters[ibody].symmetry_pars;

        // Apply constraint-based transformation
        auto& transformer = rigidbody.transformer;
        auto& param_gen = rigidbody.parameter_generator;
        auto constraint = rigidbody.constraints->get_body_constraints(ibody).at(0);

        auto new_params = param_gen->next(ibody);
        transformer->apply(std::move(new_params), constraint);

        // Verify symmetry count is preserved
        INFO("Symmetry count should be preserved after RigidTransform");
        REQUIRE(rigidbody.molecule.get_body(ibody).size_symmetry() == original_size);

        // Verify OptimizableSymmetryStorage is still there
        auto* storage = dynamic_cast<symmetry::OptimizableSymmetryStorage*>(
            rigidbody.molecule.get_body(ibody).symmetry().get_obj()
        );
        REQUIRE(storage != nullptr);

        // Undo and verify restoration
        transformer->undo();

        REQUIRE(rigidbody.molecule.get_body(ibody).size_symmetry() == original_size);
        auto& restored_sym_pars = rigidbody.conformation->absolute_parameters.parameters[ibody].symmetry_pars;
        REQUIRE(restored_sym_pars.size() == original_sym_pars.size());
    }
}

TEST_CASE("SymmetryBackup: Multiple transformations maintain symmetry integrity") {
    settings::general::verbose = false;
    settings::grid::min_bins = 250;
    settings::molecule::implicit_hydrogens = false;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::Linear;
    settings::rigidbody::iterations = 10;

    auto bodies = BodySplitter::split("tests/files/LAR1-2.pdb", {9, 99});
    bodies.get_body(0).symmetry().add(symmetry::type::p2);
    
    Rigidbody rigidbody(std::move(bodies));
    rigidbody.molecule.generate_new_hydration();
    unsigned int ibody = 0;

    auto& transformer = rigidbody.transformer;
    auto& param_gen = rigidbody.parameter_generator;

    // Apply multiple transformations in sequence
    for (int i = 0; i < 5; ++i) {
        auto params = param_gen->next(ibody);
        transformer->apply(std::move(params), ibody);

        INFO("After transformation " << i);
        REQUIRE(rigidbody.molecule.get_body(ibody).size_symmetry() == 1);
        REQUIRE(rigidbody.conformation->absolute_parameters.parameters[ibody].symmetry_pars.size() == 1);
        REQUIRE(rigidbody.conformation->initial_conformation[ibody].size_symmetry() == 1);

        // Verify OptimizableSymmetryStorage is maintained
        auto* storage = dynamic_cast<symmetry::OptimizableSymmetryStorage*>(
            rigidbody.molecule.get_body(ibody).symmetry().get_obj()
        );
        REQUIRE(storage != nullptr);
    }
}

TEST_CASE("SymmetryBackup: Grid properly sized for symmetry optimization") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::Linear;

    // Create rigidbody with symmetry
    auto bodies = BodySplitter::split("tests/files/LAR1-2.pdb", {9, 99});
    bodies.get_body(0).symmetry().add(symmetry::type::p2);
    
    Rigidbody rigidbody(std::move(bodies));
    rigidbody.molecule.generate_new_hydration();
    unsigned int ibody = 0;

    // Verify grid is large enough for initial symmetry configuration
    auto grid = rigidbody.molecule.get_grid();
    REQUIRE(grid != nullptr);
    
    // Apply a symmetry-only transformation
    auto& transformer = rigidbody.transformer;
    auto& param_gen = rigidbody.parameter_generator;
    
    // Generate parameters that modify only symmetry (no translation/rotation)
    auto params = param_gen->next(ibody);
    if (params.symmetry_pars.has_value() && !params.symmetry_pars.value().empty()) {
        // Modify symmetry parameters
        params.translation = {0, 0, 0};
        params.rotation = {0, 0, 0};
        
        // Apply the transformation
        INFO("Grid should be automatically resized to accommodate symmetry transformations");
        REQUIRE_NOTHROW(transformer->apply(std::move(params), ibody));
        
        // Verify molecule is still valid and grid contains all atoms
        REQUIRE(rigidbody.molecule.get_body(ibody).size_symmetry() == 1);
        
        // Try to regenerate hydration (this will fail if grid is too small)
        REQUIRE_NOTHROW(rigidbody.molecule.generate_new_hydration());
    }
}
