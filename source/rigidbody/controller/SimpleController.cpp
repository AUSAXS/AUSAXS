// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/controller/SimpleController.h>
#include <rigidbody/parameters/ParameterGenerationStrategy.h>
#include <rigidbody/transform/TransformStrategy.h>
#include <rigidbody/selection/BodySelectStrategy.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/constraints/ConstrainedFitter.h>
#include <rigidbody/detail/SystemSpecification.h>
#include <rigidbody/detail/MoleculeTransformParametersAbsolute.h>
#include <rigidbody/Rigidbody.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogram.h>
#include <data/Molecule.h>
#include <grid/Grid.h>

using namespace ausaxs;
using namespace ausaxs::rigidbody;
using namespace ausaxs::rigidbody::parameter;
using namespace ausaxs::rigidbody::constraints;
using namespace ausaxs::rigidbody::controller;

SimpleController::~SimpleController() = default;

void SimpleController::setup(const io::ExistingFile& measurement_path) {
    if (!calibration) {
        this->fitter = std::make_unique<fitter::ConstrainedFitter>(
            rigidbody->constraints.get(), 
            measurement_path, 
            rigidbody->molecule.get_histogram()
        );
    } else {
        auto histogram = rigidbody->molecule.get_histogram();
        histogram->apply_water_scaling_factor(calibration->get_parameter("c"));
        this->fitter = std::make_unique<fitter::ConstrainedFitter>(
            rigidbody->constraints.get(), 
            measurement_path, 
            std::move(histogram)
        );
    }
}

void SimpleController::update_fitter() {
    assert(fitter.get() != nullptr && "RigidBody::update_fitter: Fitter not initialized.");
    if (!calibration) {
        fitter->set_model(rigidbody->molecule.get_histogram());
    } else {
        auto histogram = rigidbody->molecule.get_histogram();
        histogram->apply_water_scaling_factor(calibration->get_parameter("c"));
        fitter->set_model(std::move(histogram));
    }
}

bool SimpleController::prepare_step() {
    auto& molecule = rigidbody->molecule;

    // select a body to be modified this iteration
    auto [ibody, iconstraint] = rigidbody->body_selector->next();
    if (iconstraint == -1) {    // transform free body
        auto param = rigidbody->parameter_generator->next(ibody);
        rigidbody->transformer->apply(std::move(param), ibody);
    } else {                    // transform constrained body
        DistanceConstraint& constraint = rigidbody->constraints->distance_constraints_map.at(ibody).at(iconstraint).get();
        auto param = rigidbody->parameter_generator->next(ibody);
        rigidbody->transformer->apply(std::move(param), constraint);
    }
    molecule.generate_new_hydration();

    // update the body location in the fitter
    update_fitter();
    rigidbody->conformation->absolute_parameters.chi2 = fitter->fit_chi2_only();

    step_accepted = rigidbody->conformation->absolute_parameters.chi2 < current_best_config->chi2;
    return step_accepted;
}

void SimpleController::finish_step() {
    if (step_accepted) {
        // accept the changes - update the best configuration
        *current_best_config = rigidbody->conformation->absolute_parameters;
        step_accepted = false;
    } else {
        // undo the body transforms (restores rigidbody->conformation->absolute_parameters from backup)
        rigidbody->transformer->undo();

        // regenerate grid & hydration layer
        //? potential inconsistency: the new grid may not be exactly identical to the old one,
        //? as the order of adding bodies may differ. perhaps the difference is negligible?
        rigidbody->molecule.clear_grid();
        rigidbody->molecule.generate_new_hydration();
    }
}