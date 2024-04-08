/*
This software is distributed under the GNU General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <rigidbody/RigidBody.h>
#include <rigidbody/detail/BestConf.h>
#include <rigidbody/transform/TransformFactory.h>
#include <rigidbody/selection/BodySelectFactory.h>
#include <rigidbody/parameters/ParameterGenerationFactory.h>
#include <rigidbody/constraints/ConstrainedFitter.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/parameters/Parameters.h>
#include <mini/detail/FittedParameter.h>
#include <mini/detail/Evaluation.h>
#include <utility/Exceptions.h>
#include <utility/Console.h>
#include <io/XYZWriter.h>
#include <fitter/HydrationFitter.h>
#include <fitter/LinearFitter.h>
#include <fitter/Fit.h>
#include <hydrate/Grid.h>
#include <hydrate/GridMember.h>
#include <hydrate/placement/PlacementStrategy.h>
#include <hydrate/culling/CullingStrategy.h>
#include <data/record/Atom.h>
#include <data/record/Water.h>
#include <data/Body.h>
#include <hist/intensity_calculator/DistanceHistogram.h>
#include <hist/intensity_calculator/CompositeDistanceHistogram.h>
#include <settings/RigidBodySettings.h>
#include <settings/GeneralSettings.h>
#include <plots/PlotIntensityFit.h>
#include <plots/PlotDistance.h>

using namespace rigidbody;
using namespace rigidbody::constraints;
using namespace rigidbody::parameter;

RigidBody::~RigidBody() = default;

RigidBody::RigidBody(data::Molecule&& protein) : data::Molecule(std::move(protein)) {
    initialize();
}

RigidBody::RigidBody(const data::Molecule& protein) : data::Molecule(protein) {
    initialize();
}

void RigidBody::initialize() {
    parameter_generator = factory::create_parameter_strategy(settings::rigidbody::iterations, 5, constants::pi/3);
    body_selector = factory::create_selection_strategy(this);
    transform = factory::create_transform_strategy(this);
    constraints = std::make_shared<ConstraintManager>(this);
}

void RigidBody::set_constraint_manager(std::shared_ptr<rigidbody::constraints::ConstraintManager> constraints) {
    this->constraints = constraints;
}

void RigidBody::set_body_select_manager(std::shared_ptr<rigidbody::selection::BodySelectStrategy> body_selector) {
    this->body_selector = body_selector;
}

void RigidBody::set_transform_manager(std::shared_ptr<rigidbody::transform::TransformStrategy> transform) {
    this->transform = transform;
}

void RigidBody::set_parameter_manager(std::shared_ptr<rigidbody::parameter::ParameterGenerationStrategy> parameters) {
    this->parameter_generator = parameters;
}

std::shared_ptr<fitter::Fit> RigidBody::optimize(const io::ExistingFile& measurement_path) {
    generate_new_hydration();
    prepare_fitter(measurement_path);

    if (settings::general::supplementary_plots) {
        // plots::PlotDistance::quick_plot(get_histogram(), settings::general::output + "/hist/distance_0.png");
        plots::PlotIntensityFit::quick_plot(fitter->fit().get(), settings::general::output + "initial_curve.png");
    }

    // save the best configuration in a simple struct
    detail::BestConf best(std::make_shared<grid::Grid>(*get_grid()), get_waters(), fitter->fit_chi2_only());

    if (settings::general::verbose) {
        console::print_info("\nStarting rigid body optimization.");
        std::cout << "\tInitial chi2: " << best.chi2 << std::endl;
    }

    // prepare the trajectory output
    io::XYZWriter trajectory(settings::general::output + "trajectory.xyz");
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
    update_fitter(fitter);
    auto fit = fitter->fit();
    if (calibration != nullptr) {fit->add_parameter(calibration->get_parameter("c"));}
    return fit;
}

#include <rigidbody/sequencer/All.h>
#include <rigidbody/sequencer/detail/SequenceParser.h>
std::shared_ptr<fitter::Fit> RigidBody::optimize_sequence(const io::ExistingFile& measurement_path) {
    return sequencer::SequenceParser().parse("test.txt", measurement_path, this)->execute();

    // return sequencer::Sequencer(measurement_path, this)
    //     .parameter_strategy(rigidbody::factory::create_parameter_strategy(
    //         settings::rigidbody::iterations, 
    //         5, 
    //         constants::pi/3, 
    //         settings::rigidbody::ParameterGenerationStrategyChoice::RotationsOnly
    //     ))
    //     .loop(100)
    //         .optimize()
    //             .save_on_improvement(settings::general::output + "trajectory.xyz")
    //     .end()
    //     .parameter_strategy(rigidbody::factory::create_parameter_strategy(
    //         200, 
    //         5, 
    //         constants::pi/3, 
    //         settings::rigidbody::ParameterGenerationStrategyChoice::Simple
    //     ))
    //     .loop(settings::rigidbody::iterations)
    //         .optimize()
    //             .save_on_improvement(settings::general::output + "trajectory.xyz")
    //         .end()
    //     .end()
    // .execute();
}

bool RigidBody::optimize_step(detail::BestConf& best) {
    auto grid = get_grid();

    // select a body to be modified this iteration
    auto [ibody, iconstraint] = body_selector->next();
    DistanceConstraint& constraint = constraints->distance_constraints_map.at(ibody).at(iconstraint).get();
    Parameter param = parameter_generator->next();

    Matrix R = matrix::rotation_matrix(param.alpha, param.beta, param.gamma);
    transform->apply(R, param.dr, constraint);
    generate_new_hydration(); 

    // update the body location in the fitter
    update_fitter(fitter);
    double new_chi2 = fitter->fit_chi2_only();

    // if the old configuration was better
    if (new_chi2 >= best.chi2) {
        transform->undo();          // undo the body transforms
        *grid = *best.grid;         // restore the old grid
        get_waters() = best.waters;     // restore the old waters
        return false;
    } else {
        // accept the changes
        best.grid = std::make_shared<grid::Grid>(*grid);
        best.waters = get_waters();
        best.chi2 = new_chi2;
        return true;
    }
}

void RigidBody::apply_calibration(std::shared_ptr<fitter::Fit> calibration) {
    if (settings::general::verbose) {std::cout << "\tApplying calibration to rigid body." << std::endl;}
    this->calibration = calibration;
}

void RigidBody::prepare_fitter(const std::string& measurement_path) {
    // constraints = std::make_shared<ConstraintManager>(this);
    if (calibration == nullptr) {
        fitter::ConstrainedFitter<fitter::HydrationFitter> fitter(measurement_path, get_histogram());
        fitter.set_constraint_manager(constraints);
        this->fitter = std::make_unique<fitter::ConstrainedFitter<fitter::HydrationFitter>>(std::move(fitter));
    } else {
        auto histogram = get_histogram();
        histogram->apply_water_scaling_factor(calibration->get_parameter("c"));
        fitter::ConstrainedFitter<fitter::LinearFitter> fitter(measurement_path, std::move(histogram));
        fitter.set_constraint_manager(constraints);
        this->fitter = std::make_unique<fitter::ConstrainedFitter<fitter::LinearFitter>>(std::move(fitter));
    }
}

void RigidBody::update_fitter(std::shared_ptr<fitter::LinearFitter> fitter) {
    if (calibration == nullptr) {
        fitter->set_scattering_hist(get_histogram());
    } else {
        auto histogram = get_histogram();
        histogram->apply_water_scaling_factor(calibration->get_parameter("c"));
        fitter->set_scattering_hist(std::move(histogram));
    }
}

std::shared_ptr<ConstraintManager> RigidBody::get_constraint_manager() const {
    #ifdef DEBUG
        if (constraints == nullptr) {throw except::nullptr_error("RigidBody::get_constraints: Constraint manager not initialized.");}
    #endif

    return constraints;
}