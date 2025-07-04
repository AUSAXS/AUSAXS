// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/RigidBody.h>
#include <rigidbody/detail/BestConf.h>
#include <rigidbody/transform/TransformFactory.h>
#include <rigidbody/selection/BodySelectFactory.h>
#include <rigidbody/parameters/ParameterGenerationFactory.h>
#include <rigidbody/constraints/ConstrainedFitter.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <mini/detail/FittedParameter.h>
#include <mini/detail/Evaluation.h>
#include <utility/Exceptions.h>
#include <utility/Console.h>
#include <io/detail/XYZWriter.h>
#include <fitter/SmartFitter.h>
#include <fitter/LinearFitter.h>
#include <grid/Grid.h>
#include <grid/detail/GridMember.h>
#include <data/atoms/AtomFF.h>
#include <data/atoms/Water.h>
#include <data/Body.h>
#include <hist/intensity_calculator/DistanceHistogram.h>
#include <hist/intensity_calculator/CompositeDistanceHistogram.h>
#include <settings/RigidBodySettings.h>
#include <settings/GeneralSettings.h>
#include <plots/PlotDistance.h>

using namespace ausaxs;
using namespace ausaxs::rigidbody;
using namespace ausaxs::rigidbody::constraints;
using namespace ausaxs::rigidbody::parameter;

RigidBody::~RigidBody() = default;

void RigidBody::initialize() {
    parameter_generator = factory::create_parameter_strategy(this, settings::rigidbody::iterations, 5, std::numbers::pi/3);
    body_selector = factory::create_selection_strategy(this);
    transform = factory::create_transform_strategy(this);
    constraints = std::make_shared<ConstraintManager>(this);
}

void RigidBody::set_constraint_manager(std::shared_ptr<rigidbody::constraints::ConstraintManager> constraints) {
    this->constraints = std::move(constraints);
}

void RigidBody::set_body_select_manager(std::shared_ptr<rigidbody::selection::BodySelectStrategy> body_selector) {
    this->body_selector = std::move(body_selector);
}

void RigidBody::set_transform_manager(std::shared_ptr<rigidbody::transform::TransformStrategy> transform) {
    this->transform = std::move(transform);
}

void RigidBody::set_parameter_manager(std::shared_ptr<rigidbody::parameter::ParameterGenerationStrategy> parameters) {
    this->parameter_generator = std::move(parameters);
}

std::shared_ptr<fitter::FitResult> RigidBody::optimize(const io::ExistingFile& measurement_path) {
    generate_new_hydration();
    prepare_fitter(measurement_path);

    if (settings::general::supplementary_plots) {
        // plots::PlotDistance::quick_plot(get_histogram(), settings::general::output + "/hist/distance_0.png");
        // plots::PlotIntensityFit::quick_plot(fitter->fit().get(), settings::general::output + "initial_curve.png");
    }

    // save the best configuration in a simple struct
    detail::BestConf best(std::make_shared<grid::Grid>(*get_grid()), get_waters(), fitter->fit_chi2_only());

    if (settings::general::verbose) {
        console::print_info("\nStarting rigid body optimization.");
        std::cout << "\tInitial chi2: " << best.chi2 << std::endl;
    }

    // prepare the trajectory output
    io::detail::xyz::XYZWriter trajectory(settings::general::output + "trajectory.xyz");
    trajectory.write_frame(this);

    for (unsigned int i = 0; i < settings::rigidbody::iterations; i++) {
        if (optimize_step(best)) [[unlikely]] {
            trajectory.write_frame(this);
            std::cout << "Iteration " << i << std::endl;
            console::print_success("\tRigidBody::optimize: Accepted changes. New best chi2: " + std::to_string(best.chi2));
        } else [[likely]] {
            if (i % 10 == 0 && settings::general::verbose) {
                std::cout << "Iteration " << i << "          " << std::flush;
            }
        }
    }

    save(settings::general::output + "optimized.pdb");
    update_fitter();
    auto fit = fitter->fit();
    if (calibration != nullptr) {fit->add_parameter(calibration->get_parameter("c"));}
    return fit;
}

bool RigidBody::optimize_step(detail::BestConf& best) {
    auto grid = get_grid();

    // select a body to be modified this iteration
    auto [ibody, iconstraint] = body_selector->next();
    if (iconstraint == -1) {    // transform free body
        Parameter param = parameter_generator->next(ibody);
        transform->apply(std::move(param), ibody);
    } else {                    // transform constrained body
        DistanceConstraint& constraint = constraints->distance_constraints_map.at(ibody).at(iconstraint).get();
        Parameter param = parameter_generator->next(ibody);
        transform->apply(std::move(param), constraint);
    }
    generate_new_hydration(); 

    // update the body location in the fitter
    update_fitter();
    double new_chi2 = fitter->fit_chi2_only();

    // if the old configuration was better
    if (new_chi2 >= best.chi2) {
        transform->undo();          // undo the body transforms
        *grid = *best.grid;         // restore the old grid
        get_waters() = best.waters; // restore the old waters
        signal_modified_hydration_layer();
        return false;
    } else {
        // accept the changes
        best.grid = std::make_shared<grid::Grid>(*grid);
        best.waters = get_waters();
        best.chi2 = new_chi2;
        return true;
    }
}

void RigidBody::apply_calibration(std::unique_ptr<fitter::FitResult> calibration) {
    if (settings::general::verbose) {std::cout << "\tApplying calibration to rigid body." << std::endl;}
    this->calibration = std::move(calibration);
}

void RigidBody::prepare_fitter(const io::ExistingFile& measurement_path) {
    // constraints = std::make_shared<ConstraintManager>(this);
    if (calibration == nullptr) {
        fitter::ConstrainedFitter<fitter::SmartFitter> fitter({measurement_path}, get_histogram());
        fitter.set_constraint_manager(constraints);
        this->fitter = std::make_unique<fitter::ConstrainedFitter<fitter::SmartFitter>>(std::move(fitter));
    } else {
        auto histogram = get_histogram();
        histogram->apply_water_scaling_factor(calibration->get_parameter("c"));
        fitter::ConstrainedFitter<fitter::SmartFitter> fitter(measurement_path, std::move(histogram));
        fitter.set_constraint_manager(constraints);
        this->fitter = std::make_unique<fitter::ConstrainedFitter<fitter::SmartFitter>>(std::move(fitter));
    }
}

std::unique_ptr<fitter::SmartFitter> RigidBody::get_unconstrained_fitter(const io::ExistingFile& saxs) const {
    if (calibration == nullptr) {
        return std::make_unique<fitter::SmartFitter>(saxs, get_histogram());
    } else {
        auto histogram = get_histogram();
        histogram->apply_water_scaling_factor(calibration->get_parameter("c"));
        return std::make_unique<fitter::SmartFitter>(saxs, std::move(histogram));
    }
}

std::shared_ptr<fitter::SmartFitter> RigidBody::get_fitter() const {
    return fitter;
}

void RigidBody::update_fitter() {
    if (calibration == nullptr) {
        fitter->set_model(get_histogram());
    } else {
        auto histogram = get_histogram();
        histogram->apply_water_scaling_factor(calibration->get_parameter("c"));
        fitter->set_model(std::move(histogram));
    }
}

std::shared_ptr<ConstraintManager> RigidBody::get_constraint_manager() const {
    #ifdef DEBUG
        if (constraints == nullptr) {throw except::nullptr_error("RigidBody::get_constraints: Constraint manager not initialized.");}
    #endif

    return constraints;
}