// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/controller/MetropolisController.h>
#include <rigidbody/constraints/ConstrainedFitter.h>
#include <rigidbody/parameters/ParameterGenerationStrategy.h>
#include <rigidbody/transform/TransformStrategy.h>
#include <rigidbody/selection/BodySelectStrategy.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/detail/Configuration.h>
#include <rigidbody/Rigidbody.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogram.h>
#include <data/Molecule.h>
#include <grid/Grid.h>

using namespace ausaxs;
using namespace ausaxs::rigidbody;
using namespace ausaxs::rigidbody::parameter;
using namespace ausaxs::rigidbody::constraints;
using namespace ausaxs::rigidbody::controller;

MetropolisController::~MetropolisController() = default;

void MetropolisController::setup(const io::ExistingFile& measurement_path) {
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

    auto eval_chi2 = [this] () {
        auto [ibody, iconstraint] = rigidbody->body_selector->next();
        if (iconstraint == -1) {    // transform free body
            Parameter param = rigidbody->parameter_generator->next(ibody);
            rigidbody->transformer->apply(std::move(param), ibody);
        } else {                    // transform constrained body
            DistanceConstraint& constraint = rigidbody->constraints->distance_constraints_map.at(ibody).at(iconstraint).get();
            Parameter param = rigidbody->parameter_generator->next(ibody);
            rigidbody->transformer->apply(std::move(param), constraint);
        }
        rigidbody->molecule.generate_new_hydration(); 

        // update the body location in the fitter
        update_fitter();
        return fitter->fit_chi2_only();
    };

    // initialize the best configuration with the current conformation
    rigidbody->conformation->absolute_parameters.chi2 = fitter->fit_chi2_only();
    current_best_config->chi2 = rigidbody->conformation->absolute_parameters.chi2;
}

void MetropolisController::update_fitter() {
    assert(fitter.get() != nullptr && "RigidBody::update_fitter: Fitter not initialized.");
    if (!calibration) {
        fitter->set_model(rigidbody->molecule.get_histogram());
    } else {
        auto histogram = rigidbody->molecule.get_histogram();
        histogram->apply_water_scaling_factor(calibration->get_parameter("c"));
        fitter->set_model(std::move(histogram));
    }
}

bool MetropolisController::prepare_step() {
    auto& molecule = rigidbody->molecule;

    // select a body to be modified this iteration
    auto [ibody, iconstraint] = rigidbody->body_selector->next();
    if (iconstraint == -1) {    // transform free body
        auto param = rigidbody->parameter_generator->next(ibody);
        rigidbody->transformer->apply(std::move(param), ibody);
    } else {                    // transform constrained body
        assert(iconstraint < static_cast<int>(rigidbody->constraints->get_body_constraints(ibody).size()));
        auto constraint = rigidbody->constraints->get_body_constraints(ibody)[iconstraint];
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

void MetropolisController::finish_step() {
    if (step_accepted) {
        // accept the changes - update the best configuration
        *current_best_config = rigidbody->conformation->absolute_parameters;
        current_best_config->waters = rigidbody->molecule.get_waters();
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