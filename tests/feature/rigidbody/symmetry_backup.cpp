#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <rigidbody/Rigidbody.h>
#include <rigidbody/BodySplitter.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/detail/SystemSpecification.h>
#include <rigidbody/detail/MoleculeTransformParametersAbsolute.h>
#include <rigidbody/parameters/BodyTransformParametersAbsolute.h>
#include <rigidbody/parameters/ParameterGenerationStrategy.h>
#include <rigidbody/parameters/UniformParameterGenerator.h>
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

    // Verify configuration was properly initialized with symmetry parameters
    REQUIRE(rigidbody.conformation->absolute_parameters.parameters[ibody].symmetry_pars.size() == 1);

    // Apply a transformation
    auto& transformer = rigidbody.transformer;
    auto& param_gen = rigidbody.parameter_generator;
    auto new_params = param_gen->next(ibody);
    transformer->apply(std::move(new_params), ibody);

    REQUIRE(rigidbody.molecule.get_body(ibody).size_symmetry() == 1);

    // Undo the transformation
    transformer->undo();

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
        transformer->apply(std::move(new_params), constraint, ibody);

        // Verify symmetry count is preserved
        INFO("Symmetry count should be preserved after constraint transformation");
        REQUIRE(rigidbody.molecule.get_body(ibody).size_symmetry() == original_size);

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
        transformer->apply(std::move(new_params), constraint, ibody);

        // Verify symmetry count is preserved
        INFO("Symmetry count should be preserved after RigidTransform");
        REQUIRE(rigidbody.molecule.get_body(ibody).size_symmetry() == original_size);

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

    // generate a symmetry-only perturbation
    rigidbody::parameter::UniformParameterGenerator gen(
        &rb, settings::rigidbody::iterations, {.symmetry_translation = 5, .symmetry_rotation = 0.5}
    );
    auto params = gen.next(ibody);
    REQUIRE(params.symmetry_pars.has_value());
    REQUIRE(params.symmetry_pars.value().size() == 1);

    // the generated delta must be a composite whose inner AND outer parts both received a
    // non-zero perturbation: this proves the recursion in UniformParameterGenerator
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

TEST_CASE("SymmetryBackup: undo restores the body's symmetries, not just its parameters") {
    settings::general::verbose = false;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::None;

    Molecule m{std::vector<Body>{Body(std::vector{
        AtomFF({0, 0, 0}, form_factor::form_factor_t::C), AtomFF({1, 1, 1}, form_factor::form_factor_t::C)
    })}};
    m.get_body(0).symmetry().add(symmetry::type::c2);
    Rigidbody rb(std::move(m));

    rb.parameter_generator = std::make_shared<rigidbody::parameter::UniformParameterGenerator>(
        &rb, 100, rigidbody::parameter::ParameterAmplitudes{.symmetry_rotation = 0.5}
    );

    auto body_axis = [&]() {
        return static_cast<symmetry::CyclicSymmetry*>(rb.molecule.get_body(0).symmetry().get(0))->_repeat_relation.axis;
    };
    auto param_axis = [&]() {
        auto& pars = rb.conformation->absolute_parameters.parameters[0].symmetry_pars;
        return static_cast<symmetry::CyclicSymmetry*>(pars[0].get())->_repeat_relation.axis;
    };
    auto original = body_axis();

    // a symmetry-only step takes no branch that rebuilds the body, so undo must put back the symmetries that
    // apply_symmetry_delta wrote into it. Discarding the body backup there would leave it holding the rejected axis
    // while the parameters held the restored one, and every later delta would be computed against the wrong state.
    for (int i = 0; i < 20; ++i) {
        rb.transformer->apply(rb.parameter_generator->next(0), 0);
        rb.transformer->undo();

        REQUIRE_THAT((body_axis() - param_axis()).magnitude(), Catch::Matchers::WithinAbs(0, 1e-9));
        REQUIRE_THAT((body_axis() - original).magnitude(), Catch::Matchers::WithinAbs(0, 1e-9));
    }
}

TEST_CASE("SymmetryBackup: undo restores the body's symmetries under constrained transforms too") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::grid::min_bins = 250;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::Backbone;

    auto strategy = GENERATE(
        settings::rigidbody::TransformationStrategyChoice::SingleTransform,
        settings::rigidbody::TransformationStrategyChoice::RigidTransform
    );
    settings::rigidbody::transform_strategy = strategy;

    auto bodies = BodySplitter::split("tests/files/LAR1-2.pdb", {9, 99});
    bodies.get_body(0).symmetry().add(symmetry::type::c2);
    Rigidbody rb(std::move(bodies));
    rb.molecule.generate_new_hydration();
    unsigned int ibody = 0;

    // symmetry-only steps take the branch that leaves the bodies untouched, which is exactly where the discarded
    // backup used to strand them: the parameters were rolled back while the body kept the rejected axis.
    rb.parameter_generator = std::make_shared<rigidbody::parameter::UniformParameterGenerator>(
        &rb, 100, rigidbody::parameter::ParameterAmplitudes{.symmetry_rotation = 0.5}
    );
    auto constraint = rb.constraints->get_body_constraints(ibody).at(0);

    auto body_axis = [&]() {
        return static_cast<symmetry::CyclicSymmetry*>(rb.molecule.get_body(ibody).symmetry().get(0))->_repeat_relation.axis;
    };
    auto param_axis = [&]() {
        auto& pars = rb.conformation->absolute_parameters.parameters[ibody].symmetry_pars;
        return static_cast<symmetry::CyclicSymmetry*>(pars[0].get())->_repeat_relation.axis;
    };
    // move the bodies off the initial conformation and keep the result. Without this the destructive reset that
    // RigidTransform used to perform before the branch would be indistinguishable from a no-op, since the bodies
    // still sit exactly where initial_conformation put them.
    {
        rigidbody::parameter::UniformParameterGenerator warmup(
            &rb, 100, rigidbody::parameter::ParameterAmplitudes{.translation = 2, .rotation = 0.5}
        );
        rb.transformer->apply(warmup.next(ibody), constraint, ibody);
    }

    auto original = body_axis();

    // RigidTransform moves a whole group of bodies, and the body the symmetry delta belongs to may sit outside that
    // group, so checking only its axis would miss a group left stranded in the initial conformation. Snapshot the
    // coordinates of every body as well: undo must be a no-op over the entire molecule.
    auto snapshot = [&]() {
        std::vector<Vector3<double>> out;
        for (unsigned int i = 0; i < rb.molecule.size_body(); ++i) {
            for (const auto& a : rb.molecule.get_body(i).get_atoms()) {out.push_back(a.coordinates());}
        }
        return out;
    };
    auto original_coords = snapshot();

    for (int i = 0; i < 10; ++i) {
        rb.transformer->apply(rb.parameter_generator->next(ibody), constraint, ibody);
        rb.transformer->undo();

        REQUIRE_THAT((body_axis() - param_axis()).magnitude(), Catch::Matchers::WithinAbs(0, 1e-9));
        REQUIRE_THAT((body_axis() - original).magnitude(), Catch::Matchers::WithinAbs(0, 1e-9));

        auto coords = snapshot();
        REQUIRE(coords.size() == original_coords.size());
        double worst = 0;
        for (std::size_t j = 0; j < coords.size(); ++j) {worst = std::max(worst, (coords[j] - original_coords[j]).magnitude());}
        REQUIRE_THAT(worst, Catch::Matchers::WithinAbs(0, 1e-9));
    }
}
