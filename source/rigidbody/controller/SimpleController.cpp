// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/controller/SimpleController.h>

#include <data/Molecule.h>
#include <dataset/SimpleDataset.h>
#include <hydrate/ExplicitHydration.h>
#include <rigidbody/Rigidbody.h>
#include <rigidbody/constraints/ConstrainedFitter.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/detail/MoleculeTransformParametersAbsolute.h>
#include <rigidbody/detail/SystemSpecification.h>
#include <rigidbody/parameters/ParameterGenerationStrategy.h>  // IWYU pragma: keep
#include <rigidbody/selection/BodySelectStrategy.h>  // IWYU pragma: keep
#include <rigidbody/selection/ParameterMask.h>
#include <rigidbody/transform/TransformStrategy.h>  // IWYU pragma: keep
#include <settings/HistogramSettings.h>

using namespace ausaxs;
using namespace ausaxs::rigidbody;
using namespace ausaxs::rigidbody::parameter;
using namespace ausaxs::rigidbody::constraints;
using namespace ausaxs::rigidbody::controller;

namespace {
    // replace the hydration shell with a snapshot from get_waters(); copied, since the snapshot must outlive further rejected steps
    void restore_hydration(data::Molecule& molecule, const std::vector<data::Water>& waters) {
        molecule.clear_hydration();
        auto hydration = std::make_unique<hydrate::ExplicitHydration>(waters);
        hydration->expanded_across_symmetry = true; // snapshotted from a molecule already hydrated by generate_new_hydration()
        molecule.get_body(0).set_hydration(std::move(hydration));
    }
}

SimpleController::~SimpleController() = default;

void SimpleController::setup(const io::ExistingFile& measurement_path) {
    SimpleDataset data(measurement_path);
    if (settings::axes::rebin) {data.rebin();}
    this->fitter = std::make_unique<fitter::ConstrainedFitter>(
        rigidbody->constraints.get(),
        std::move(data),
        rigidbody->molecule.get_histogram()
    );

    // initialize the best configuration with the current conformation
    rigidbody->conformation->absolute_parameters.chi2 = fitter->fit_chi2_only();
    *current_best_config = rigidbody->conformation->absolute_parameters;
    current_best_config->waters = rigidbody->molecule.get_waters(); // a step rejected before any is accepted restores these
}

void SimpleController::update_fitter() {
    assert(fitter.get() != nullptr && "RigidBody::update_fitter: Fitter not initialized.");
    fitter->set_model(rigidbody->molecule.get_histogram());
}

bool SimpleController::prepare_step() {
    auto& molecule = rigidbody->molecule;

    // select a body to be modified this iteration
    auto [ibody, iconstraint, mask] = rigidbody->body_selector->next_mask();
    if (iconstraint == -1) {    // transform free body
        auto param = rigidbody->parameter_generator->next(ibody);
        mask.apply(param);
        rigidbody->transformer->apply(param, ibody);
    } else {                    // transform constrained body
        assert(iconstraint < static_cast<int>(rigidbody->constraints->get_body_constraints(ibody).size()));
        auto* constraint = rigidbody->constraints->get_body_constraints(ibody)[iconstraint];
        auto param = rigidbody->parameter_generator->next(ibody);
        mask.apply(param);
        rigidbody->transformer->apply(param, constraint, ibody);
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
        current_best_config->waters = rigidbody->molecule.get_waters();
        step_accepted = false;
    } else {
        // undo the body transforms (restores rigidbody->conformation->absolute_parameters from backup)
        rigidbody->transformer->undo();
        rigidbody->conformation->absolute_parameters.chi2 = current_best_config->chi2; // also restore chi2

        // the grid still holds the rejected body positions, so it is rebuilt from the restored ones on next use
        rigidbody->molecule.clear_grid();

        // restore the best conformation's hydration shell rather than generating a new one
        if (current_best_config->waters.empty()) {
            rigidbody->molecule.generate_new_hydration(); // no shell to restore, e.g. a hydration strategy that adds no waters
        } else {
            restore_hydration(rigidbody->molecule, current_best_config->waters);
        }
    }
}