#include <rigidbody/RigidBody.h>

#include <Symbols.h>
#include <rigidbody/transform/TransformFactory.h>
#include <rigidbody/selection/BodySelectFactory.h>
#include <rigidbody/parameters/ParameterGenerationFactory.h>
#include <rigidbody/RigidBodySettings.h>
#include <utility/Exceptions.h>
#include <io/XYZWriter.h>
#include <plots/PlotIntensityFit.h>
#include <plots/PlotDistance.h>

using namespace rigidbody;

RigidBody::RigidBody(Protein&& protein) : Protein(std::move(protein)) {
    setup();
}

RigidBody::RigidBody(const Protein& protein) : Protein(protein) {
    setup();
}

void RigidBody::setup() {
    constraints = factory::create_parameter_strategy(settings::rigidbody::iterations, 5, M_PI/3, settings::rigidbody::culling_strategy);
}

std::shared_ptr<fitter::Fit> RigidBody::optimize(std::string measurement_path) {
    generate_new_hydration();
    auto fitter = prepare_fitter(measurement_path);
    double best_chi2 = fitter->fit_only();
    plots::PlotDistance::quick_plot(get_histogram(), settings::general::output + "/hist/distance_0.png");

    if (settings::general::supplementary_plots) {
        plots::PlotIntensityFit::quick_plot(fitter->fit(), settings::general::output + "initial_curve.png");
    }

    if (settings::general::verbose) {
        utility::print_info("\nStarting rigid body optimization.");
        std::cout << "\tInitial chi2: " << best_chi2 << std::endl;
    }

    Parameters params(this);
    io::XYZWriter writer(settings::general::output + "trajectory.xyz");
    writer.write_frame(this);

    // save the best configuration so we can restore it after each failed attempt
    std::shared_ptr<Grid> grid = get_grid();
    std::shared_ptr<Grid> best_grid = std::make_shared<Grid>(*grid);
    std::vector<Water> best_waters = waters();
    unsigned int optimized_step = 0;
    for (unsigned int i = 0; i < settings::rigidbody::iterations; i++) {
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
        if (new_chi2 >= best_chi2) {
            transform->undo();          // undo the body transforms
            *grid = *best_grid;         // restore the old grid
            waters() = best_waters;     // restore the old waters

            if (i % 10 == 0 && settings::general::verbose) {
                std::cout << "\rIteration " << i << "          " << std::flush;
            }
        } else {
            writer.write_frame(this);

            // accept the changes
            best_grid = std::make_shared<Grid>(*grid);
            best_waters = waters();
            best_chi2 = new_chi2;
            // params.update(body.uid, param);
            std::cout << "\rIteration " << i << std::endl;
            utility::print_success("\tRigidBody::optimize: Accepted changes. New best chi2: " + std::to_string(new_chi2));
            // plots::PlotDistance::quick_plot(get_histogram(), settings::general::output + "/hist/distance_" + std::to_string(++optimized_step) + ".png");
        }
    }

    save(settings::general::output + "optimized.pdb");
    update_fitter(fitter);
    auto fit = fitter->fit();
    if (calibration != nullptr) {fit->add_parameter(calibration->get_parameter("c"));}
    return fit;
}

void RigidBody::apply_calibration(std::shared_ptr<fitter::Fit> calibration) {
    if (settings::general::verbose) {std::cout << "\tApplying calibration to rigid body." << std::endl;}
    this->calibration = calibration;
}

std::shared_ptr<fitter::LinearFitter> RigidBody::prepare_fitter(std::string measurement_path) {
    constraints = std::make_shared<ConstraintManager>(this);
    if (calibration == nullptr) {
        fitter::ConstrainedFitter<fitter::HydrationFitter> fitter(measurement_path, get_histogram());
        fitter.set_constraint_manager(constraints);
        return std::make_shared<fitter::ConstrainedFitter<fitter::HydrationFitter>>(std::move(fitter));
    } else {
        auto histogram = get_histogram();
        histogram.apply_water_scaling_factor(calibration->get_parameter("c"));
        fitter::ConstrainedFitter<fitter::LinearFitter> fitter(measurement_path, std::move(histogram));
        fitter.set_constraint_manager(constraints);
        return std::make_shared<fitter::ConstrainedFitter<fitter::LinearFitter>>(std::move(fitter));
    }
}

void RigidBody:: update_fitter(std::shared_ptr<fitter::LinearFitter> fitter) {
    if (calibration == nullptr) {
        fitter->set_scattering_hist(get_histogram());
    } else {
        auto histogram = get_histogram();
        histogram.apply_water_scaling_factor(calibration->get_parameter("c"));
        fitter->set_scattering_hist(std::move(histogram));
    }
}
