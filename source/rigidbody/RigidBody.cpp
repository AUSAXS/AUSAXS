#include <rigidbody/RigidBody.h>

#include <Symbols.h>
#include <rigidbody/transform/TransformFactory.h>
#include <rigidbody/selection/BodySelectFactory.h>
#include <rigidbody/parameters/ParameterGenerationFactory.h>
#include <rigidbody/constraints/ConstrainedFitter.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/parameters/Parameters.h>
#include <mini/detail/FittedParameter.h>
#include <utility/Exceptions.h>
#include <utility/Console.h>
#include <io/XYZWriter.h>
#include <plots/PlotIntensityFit.h>
#include <plots/PlotDistance.h>
#include <settings/RigidBodySettings.h>
#include <settings/GeneralSettings.h>
#include <fitter/HydrationFitter.h>
#include <fitter/LinearFitter.h>
#include <fitter/Fit.h>

using namespace rigidbody;

RigidBody::~RigidBody() = default;

RigidBody::RigidBody(Protein&& protein) : Protein(std::move(protein)) {
    initialize();
}

RigidBody::RigidBody(const Protein& protein) : Protein(protein) {
    initialize();
}

void RigidBody::initialize() {
    parameter_generator = factory::create_parameter_strategy(settings::rigidbody::iterations, 5, M_PI/3);
    body_selector = factory::create_selection_strategy(this);
    transform = factory::create_transform_strategy(this);
}

std::shared_ptr<fitter::Fit> RigidBody::optimize(const std::string& measurement_path) {
    generate_new_hydration();
    prepare_fitter(measurement_path);

    if (settings::general::supplementary_plots) {
        // plots::PlotDistance::quick_plot(get_histogram(), settings::general::output + "/hist/distance_0.png");
        plots::PlotIntensityFit::quick_plot(fitter->fit(), settings::general::output + "initial_curve.png");
    }

    // save the best configuration in a simple struct
    detail::BestConf best {
        .grid = std::make_shared<grid::Grid>(*get_grid()),
        .waters = waters(),
        .chi2 = fitter->fit_only()
    };

    if (settings::general::verbose) {
        console::print_info("\nStarting rigid body optimization.");
        std::cout << "\tInitial chi2: " << best.chi2 << std::endl;
    }

    // prepare the trajectory output
    io::XYZWriter trajectory(settings::general::output + "trajectory.xyz");
    trajectory.write_frame(this);

    unsigned int optimized_step = 0;
    for (unsigned int i = 0; i < settings::rigidbody::iterations; i++) {
        if (optimize_step(best)) [[unlikely]] {
            trajectory.write_frame(this);
            std::cout << "\rIteration " << i << std::endl;
            console::print_success("\tRigidBody::optimize: Accepted changes. New best chi2: " + std::to_string(best.chi2));
        } else [[likely]] {
            if (i % 10 == 0 && settings::general::verbose) {
                std::cout << "\rIteration " << i << "          " << std::flush;
            }
        }
    }

    save(settings::general::output + "optimized.pdb");
    update_fitter(fitter);
    auto fit = fitter->fit();
    if (calibration != nullptr) {fit->add_parameter(calibration->get_parameter("c"));}
    return fit;
}

bool RigidBody::optimize_step(detail::BestConf& best) {
    std::shared_ptr<grid::Grid> grid = get_grid();

    // select a body to be modified this iteration
    auto [ibody, iconstraint] = body_selector->next();
    std::shared_ptr<DistanceConstraint> constraint = constraints->distance_constraints_map.at(ibody).at(iconstraint);
    Parameter param = parameter_generator->next();

    Matrix R = matrix::rotation_matrix(param.alpha, param.beta, param.gamma);
    transform->apply(R, param.dx, constraint);
    generate_new_hydration(); 

    // update the body location in the fitter
    update_fitter(fitter);
    double new_chi2 = fitter->fit_only();

    // if the old configuration was better
    if (new_chi2 >= best.chi2) {
        transform->undo();          // undo the body transforms
        *grid = *best.grid;         // restore the old grid
        waters() = best.waters;     // restore the old waters
    } else {
        // accept the changes
        best.grid = std::make_shared<grid::Grid>(*grid);
        best.waters = waters();
        best.chi2 = new_chi2;
    }
}

void RigidBody::apply_calibration(std::shared_ptr<fitter::Fit> calibration) {
    if (settings::general::verbose) {std::cout << "\tApplying calibration to rigid body." << std::endl;}
    this->calibration = calibration;
}

void RigidBody::prepare_fitter(const std::string& measurement_path) {
    constraints = std::make_shared<ConstraintManager>(this);
    if (calibration == nullptr) {
        fitter::ConstrainedFitter<fitter::HydrationFitter> fitter(measurement_path, get_histogram());
        fitter.set_constraint_manager(constraints);
        this->fitter = std::make_shared<fitter::ConstrainedFitter<fitter::HydrationFitter>>(std::move(fitter));
    } else {
        auto histogram = get_histogram();
        histogram.apply_water_scaling_factor(calibration->get_parameter("c"));
        fitter::ConstrainedFitter<fitter::LinearFitter> fitter(measurement_path, std::move(histogram));
        fitter.set_constraint_manager(constraints);
        this->fitter = std::make_shared<fitter::ConstrainedFitter<fitter::LinearFitter>>(std::move(fitter));
    }
}

void RigidBody::update_fitter(std::shared_ptr<fitter::LinearFitter> fitter) {
    if (calibration == nullptr) {
        fitter->set_scattering_hist(get_histogram());
    } else {
        auto histogram = get_histogram();
        histogram.apply_water_scaling_factor(calibration->get_parameter("c"));
        fitter->set_scattering_hist(std::move(histogram));
    }
}