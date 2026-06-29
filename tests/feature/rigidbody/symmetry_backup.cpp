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
#include <rigidbody/parameters/ParameterGenerationStrategies.h>
#include <rigidbody/transform/TransformStrategy.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <data/symmetry/BodySymmetryFacade.h>
#include <data/symmetry/CyclicSymmetry.h>
#include <data/symmetry/CompositeSymmetry.h>
#include <settings/All.h>
#include <io/ExistingFile.h>

#include <algorithm>
#include <numbers>
#include <span>
#include <vector>

using namespace ausaxs;
using namespace ausaxs::data;
using namespace ausaxs::rigidbody;

TEST_CASE("SymmetryBackup: Symmetry structure preserved in original_conformation") {
    settings::general::verbose = false;
    settings::grid::min_bins = 250;
    settings::molecule::implicit_hydrogens = false;

    SECTION("via direct construction") {
        auto bodies = BodySplitter::split("tests/files/LAR1-2.pdb", {9, 99});
        bodies.get_body(0).symmetry().add(symmetry::type::c2);
        
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
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::Backbone;

    // Create a rigidbody with symmetry
    auto bodies = BodySplitter::split("tests/files/LAR1-2.pdb", {9, 99});
    bodies.get_body(0).symmetry().add(symmetry::type::c2);
    
    Rigidbody rigidbody(std::move(bodies));
    rigidbody.molecule.generate_new_hydration();
    unsigned int ibody = 0;

    // Store original symmetry parameters
    auto original_sym_pars = rigidbody.conformation->absolute_parameters.parameters[ibody];
    REQUIRE(original_sym_pars.symmetry_pars.size() == 1);

    // Apply a transformation that modifies symmetry parameters
    auto& transformer = rigidbody.transformer;
    auto& param_gen = rigidbody.parameter_generator;

    // Generate new parameters (which should include modified symmetry parameters)
    auto new_params = param_gen->next(ibody);

    // Store the new symmetry parameters for verification
    auto expected_new_sym_pars = new_params;
    REQUIRE(expected_new_sym_pars.symmetry_pars.has_value());
    REQUIRE(expected_new_sym_pars.symmetry_pars.value().size() == 1);

    // Apply the transformation
    transformer->apply(std::move(new_params), ibody);

    // Verify symmetry parameters were updated in configuration
    auto& updated_sym_pars = rigidbody.conformation->absolute_parameters.parameters[ibody].symmetry_pars;
    REQUIRE(updated_sym_pars.size() == 1);

    INFO("Symmetry parameters should be updated after transformation");
    auto t1 = updated_sym_pars[0]->span_translation();
    auto t2 = expected_new_sym_pars.symmetry_pars.value()[0]->span_translation();
    for (int i = 0; i < 3; ++i) {REQUIRE_THAT(t1[i], Catch::Matchers::WithinAbs(t2[i], 1e-6));}

    // Undo the transformation
    transformer->undo();

    // Verify symmetry parameters were restored
    auto& restored_sym_pars = rigidbody.conformation->absolute_parameters.parameters[ibody].symmetry_pars;
    REQUIRE(restored_sym_pars.size() == 1);

    INFO("Symmetry parameters should be restored after undo");
    t1 = restored_sym_pars[0]->span_translation();
    t2 = original_sym_pars.symmetry_pars[0]->span_translation();
    for (int i = 0; i < 3; ++i) {REQUIRE_THAT(t1[i], Catch::Matchers::WithinAbs(t2[i], 1e-6));}
    t1 = restored_sym_pars[0]->span_rotation();
    t2 = original_sym_pars.symmetry_pars[0]->span_rotation();
    for (int i = 0; i < 3; ++i) {REQUIRE_THAT(t1[i], Catch::Matchers::WithinAbs(t2[i], 1e-6));} 
}

TEST_CASE("SymmetryBackup: Body symmetry storage preserved through transformations") {
    settings::general::verbose = false;
    settings::grid::min_bins = 250;
    settings::molecule::implicit_hydrogens = false;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::Backbone;

    auto bodies = BodySplitter::split("tests/files/LAR1-2.pdb", {9, 99});
    bodies.get_body(0).symmetry().add(symmetry::type::c2);
    
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
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::Backbone;

    SECTION("SingleTransform") {
        settings::grid::min_bins = 250;
        settings::rigidbody::transform_strategy = settings::rigidbody::TransformationStrategyChoice::SingleTransform;

        auto bodies = BodySplitter::split("tests/files/LAR1-2.pdb", {9, 99});
        bodies.get_body(0).symmetry().add(symmetry::type::c2);
        
        Rigidbody rigidbody(std::move(bodies));
        rigidbody.molecule.generate_new_hydration();
        unsigned int ibody = 0;

        // Store original symmetry state
        auto original_size = rigidbody.molecule.get_body(ibody).size_symmetry();
        auto original_sym_pars = rigidbody.conformation->absolute_parameters.parameters[ibody];

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
        REQUIRE(restored_sym_pars.size() == original_sym_pars.symmetry_pars.size());
    }

    SECTION("RigidTransform") {
        settings::grid::min_bins = 250;
        settings::rigidbody::transform_strategy = settings::rigidbody::TransformationStrategyChoice::RigidTransform;

        auto bodies = BodySplitter::split("tests/files/LAR1-2.pdb", {9, 99});
        bodies.get_body(0).symmetry().add(symmetry::type::c2);
        
        Rigidbody rigidbody(std::move(bodies));
        rigidbody.molecule.generate_new_hydration();
        unsigned int ibody = 0;

        // Store original symmetry state
        auto original_size = rigidbody.molecule.get_body(ibody).size_symmetry();
        auto original_sym_pars = rigidbody.conformation->absolute_parameters.parameters[ibody];

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
        auto& restored_sym_pars = rigidbody.conformation->absolute_parameters.parameters[ibody];
        REQUIRE(restored_sym_pars.symmetry_pars.size() == original_sym_pars.symmetry_pars.size());
    }
}

TEST_CASE("SymmetryBackup: Multiple transformations maintain symmetry integrity") {
    settings::general::verbose = false;
    settings::grid::min_bins = 250;
    settings::molecule::implicit_hydrogens = false;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::Backbone;
    settings::rigidbody::iterations = 10;

    auto bodies = BodySplitter::split("tests/files/LAR1-2.pdb", {9, 99});
    bodies.get_body(0).symmetry().add(symmetry::type::c2);
    
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

TEST_CASE("SymmetryBackup: CompositeSymmetry parameters are optimised") {
    settings::general::verbose = false;
    settings::grid::min_bins = 100;
    settings::molecule::implicit_hydrogens = false;
    settings::molecule::center = false;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::None;

    // a p2-3-style composite: an inner c2 (with an offset) nested inside an outer c3
    auto make_composite = [] {
        auto inner = std::make_unique<symmetry::CyclicSymmetry>(
            Vector3<double>{4, 0, 0}, Vector3<double>{0, 0, 0}, Vector3<double>{0, 0, 1}, std::numbers::pi, 1
        );
        auto outer = std::make_unique<symmetry::CyclicSymmetry>(
            Vector3<double>{0, 0, 0}, Vector3<double>{0, 0, 0}, Vector3<double>{0, 0, 1}, 2*std::numbers::pi/3, 2
        );
        return std::make_unique<symmetry::CompositeSymmetry>(std::move(inner), std::move(outer));
    };

    Molecule m({Body{std::vector{
        AtomFF({1, 0, 0}, form_factor::form_factor_t::C),
        AtomFF({0, 1, 0}, form_factor::form_factor_t::C),
        AtomFF({0, 0, 1}, form_factor::form_factor_t::C)
    }}});
    m.get_body(0).symmetry().add(make_composite());

    Rigidbody rb(std::move(m));
    unsigned int ibody = 0;
    REQUIRE(rb.molecule.get_body(ibody).size_symmetry() == 1);
    REQUIRE(rb.conformation->absolute_parameters.parameters[ibody].symmetry_pars.size() == 1);

    // enable optimisation of the symmetry parameters (composites cannot go through the
    // type-based add() that normally sets these flags, so set them directly)
    auto* storage = dynamic_cast<symmetry::OptimizableSymmetryStorage*>(rb.molecule.get_body(ibody).symmetry().get_obj());
    REQUIRE(storage != nullptr);
    storage->optimize_translate = true;
    storage->optimize_rot_axis = true;

    // generate a symmetry-only perturbation
    rigidbody::parameter::SymmetryOnly gen(&rb, settings::rigidbody::iterations, 5, 0.5);
    auto params = gen.next(ibody);
    REQUIRE(params.symmetry_pars.has_value());
    REQUIRE(params.symmetry_pars.value().size() == 1);

    // the generated delta must be a composite whose inner AND outer parts both received a
    // non-zero perturbation: this proves the recursion in ParameterGenerationStrategies
    auto* delta = dynamic_cast<symmetry::CompositeSymmetry*>(params.symmetry_pars.value()[0].get());
    REQUIRE(delta != nullptr);
    auto nonzero = [](std::span<double> s) {
        return std::any_of(s.begin(), s.end(), [](double v) {return v != 0;});
    };
    CHECK(nonzero(delta->inner->span_translation()));
    CHECK(nonzero(delta->inner->span_rotation()));
    CHECK(nonzero(delta->outer->span_translation()));
    CHECK(nonzero(delta->outer->span_rotation()));

    // capture the inner/outer offsets before applying, so we can confirm they change
    auto* live = dynamic_cast<symmetry::CompositeSymmetry*>(rb.molecule.get_body(ibody).symmetry().get(0));
    REQUIRE(live != nullptr);
    std::vector<double> inner_before(live->inner->span_translation().begin(), live->inner->span_translation().end());

    // apply: this exercises the recursion in TransformStrategy (add_symmetries + apply_symmetry)
    rb.transformer->apply(std::move(params), ibody);

    // apply move-reassigns the body, so re-fetch the live symmetry object
    live = dynamic_cast<symmetry::CompositeSymmetry*>(rb.molecule.get_body(ibody).symmetry().get(0));
    REQUIRE(live != nullptr);

    // the live symmetry's inner offset must have changed, and must equal the accumulated
    // absolute parameters that were recorded for the conformation
    auto* accumulated = dynamic_cast<symmetry::CompositeSymmetry*>(
        rb.conformation->absolute_parameters.parameters[ibody].symmetry_pars[0].get()
    );
    REQUIRE(accumulated != nullptr);

    auto live_inner_t = live->inner->span_translation();
    auto acc_inner_t = accumulated->inner->span_translation();
    bool changed = false;
    for (int i = 0; i < 3; ++i) {
        REQUIRE_THAT(live_inner_t[i], Catch::Matchers::WithinAbs(acc_inner_t[i], 1e-9));
        if (live_inner_t[i] != inner_before[i]) {changed = true;}
    }
    CHECK(changed);
}

TEST_CASE("SymmetryBackup: Grid properly sized for symmetry optimization") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::Backbone;

    // Create rigidbody with symmetry
    auto bodies = BodySplitter::split("tests/files/LAR1-2.pdb", {9, 99});
    bodies.get_body(0).symmetry().add(symmetry::type::c2);
    
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
