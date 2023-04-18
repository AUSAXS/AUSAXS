#include <Symbols.h>
#include <rigidbody/RigidBody.h>
#include <rigidbody/parameters/Parameters.h>
#include <rigidbody/parameters/SimpleParameterGeneration.h>
#include <rigidbody/parameters/RotationsOnly.h>
#include <rigidbody/transform/RigidTransform.h>
#include <rigidbody/transform/SingleTransform.h>
#include <rigidbody/selection/SequentialSelect.h>
#include <rigidbody/selection/RandomSelect.h>
#include <rigidbody/selection/RandomConstraintSelect.h>
#include <rigidbody/ConstrainedFitter.h>
#include <rigidbody/DistanceConstraint.h>
#include <rigidbody/OverlapConstraint.h>
#include <utility/Exceptions.h>
#include <math/Matrix.h>
#include <math/MatrixUtils.h>
#include <io/XYZWriter.h>
#include <plots/PlotIntensityFit.h>
#include <plots/PlotDistance.h>

#include <sstream>

using namespace rigidbody;

RigidBody::RigidBody(Protein&& protein) : Protein(std::move(protein)) {
    setup();
}

RigidBody::RigidBody(const Protein& protein) : Protein(protein) {
    setup();
}

void RigidBody::setup() {
    // Set body transformation strategy
    switch (setting::rigidbody::tsc) {
        case setting::rigidbody::TransformationStrategyChoice::RigidTransform:
            transform = std::make_unique<RigidTransform>(this); 
            break;
        case setting::rigidbody::TransformationStrategyChoice::SingleTransform:
            transform = std::make_unique<SingleTransform>(this);
            break;
        default: 
            throw except::unknown_argument("RigidBody::RigidBody: Unkown TransformationStrategy.");
    }

    // Set parameter generation strategy
    switch (setting::rigidbody::pgsc) {
        case setting::rigidbody::ParameterGenerationStrategyChoice::Simple:
            parameter_generator = std::make_unique<SimpleParameterGeneration>(setting::rigidbody::iterations, 5, M_PI/3);
            break;
        case setting::rigidbody::ParameterGenerationStrategyChoice::RotationsOnly:
            parameter_generator = std::make_unique<RotationsOnly>(setting::rigidbody::iterations, 5, M_PI/3);
            break;
        default: 
            throw except::unknown_argument("RigidBody::RigidBody: Unknown ParameterGenerationStrategy.");
    }

    // Set body selection strategy
    switch (setting::rigidbody::bssc) {
        case setting::rigidbody::BodySelectStrategyChoice::RandomSelect:
            body_selector = std::make_unique<RandomSelect>(this);
            break;
        case setting::rigidbody::BodySelectStrategyChoice::RandomConstraintSelect:
            body_selector = std::make_unique<RandomConstraintSelect>(this);
            break;
        case setting::rigidbody::BodySelectStrategyChoice::SequentialSelect:
            body_selector = std::make_unique<SequentialSelect>(this);
            break;
        default: 
            throw except::unknown_argument("RigidBody::RigidBody: Unknown BodySelectStrategy.");
    }
}

std::shared_ptr<fitter::Fit> RigidBody::optimize(std::string measurement_path) {
    generate_new_hydration();
    generate_constraint_map();
    auto fitter = prepare_fitter(measurement_path);

    double best_chi2 = fitter->fit_only();

    if (setting::general::supplementary_plots) {
        plots::PlotIntensityFit::quick_plot(fitter->fit(), setting::general::output + "initial_curve.png");
    }

    if (setting::general::verbose) {
        utility::print_info("\nStarting rigid body optimization.");
        std::cout << "\tInitial chi2: " << best_chi2 << std::endl;
    }

    Parameters params(this);
    io::XYZWriter writer(setting::general::output + "trajectory.xyz");
    writer.write_frame(this);

    // save the best configuration so we can restore it after each failed attempt
    std::shared_ptr<Grid> grid = get_grid();
    std::shared_ptr<Grid> best_grid = std::make_shared<Grid>(*grid);
    std::vector<Water> best_waters = waters();
    unsigned int optimized_step = 0;
    for (unsigned int i = 0; i < setting::rigidbody::iterations; i++) {
        // select a body to be modified this iteration
        auto [ibody, iconstraint] = body_selector->next();
        std::shared_ptr<Constraint> constraint = constraint_map.at(ibody).at(iconstraint);
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

            if (i % 10 == 0 && setting::general::verbose) {
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
            plots::PlotDistance::quick_plot(get_histogram(), setting::general::output + "/hist/distance_" + std::to_string(++optimized_step) + ".png");
        }
    }

    save(setting::general::output + "optimized.pdb");
    update_fitter(fitter);
    auto fit = fitter->fit();
    if (calibration != nullptr) {fit->add_parameter(calibration->get_parameter("c"));}
    return fit;
}

void RigidBody::generate_simple_volume_constraints() {
    if (setting::general::verbose) {utility::print_info("\tGenerating simple constraints for rigid body optimization.");}
    for (unsigned int ibody1 = 0; ibody1 < bodies.size(); ibody1++) {
        for (unsigned int ibody2 = ibody1+1; ibody2 < bodies.size(); ibody2++) {
            const Body& body1 = body(ibody1);
            const Body& body2 = body(ibody2);

            double min_dist = std::numeric_limits<double>::max();
            int min_atom1 = -1, min_atom2 = -1;
            for (unsigned int iatom1 = 0; iatom1 < body1.atoms().size(); iatom1++) {
                const Atom& atom1 = body1.atoms(iatom1);
                if (atom1.element != constants::symbols::carbon) {continue;}

                for (unsigned int iatom2 = 0; iatom2 < body2.atoms().size(); iatom2++) {
                    const Atom& atom2 = body2.atoms(iatom2);
                    if (atom2.element != constants::symbols::carbon) {continue;}

                    double dist = atom1.distance(atom2);
                    if (dist > min_dist) {continue;}

                    min_dist = dist;
                    min_atom1 = iatom1;
                    min_atom2 = iatom2;
                }
            }

            // no carbon atoms found
            if (min_atom1 == -1 || min_atom2 == -1) {continue;}

            // check if the bodies are close enough for a constraint to make sense
            if (min_dist > setting::rigidbody::bond_distance) {continue;} 
            rigidbody::DistanceConstraint constraint(this, ibody1, ibody2, min_atom1, min_atom2);
            add_constraint(std::move(constraint));

            if (setting::general::verbose) {
                std::cout << "\tConstraint created between bodies " << ibody1 << " and " << ibody2 << " on atoms " << body1.atoms(min_atom1).name << " and " << body2.atoms(min_atom2).name << std::endl;
            }
        }
    }
    if (constraints.empty()) {
        throw except::unexpected("RigidBody::generate_simple_constraints: No constraints were generated. This is probably a bug.");
    }

    generate_constraint_map();
}

void RigidBody::generate_simple_linear_constraints() {
    if (setting::general::verbose) {utility::print_info("\tGenerating simple constraints for rigid body optimization.");}
    for (unsigned int ibody1 = 0; ibody1 < bodies.size()-1; ibody1++) {
        unsigned int ibody2 = ibody1 + 1;

        const Body& body1 = body(ibody1);
        const Body& body2 = body(ibody2);

        double min_dist = std::numeric_limits<double>::max();
        int min_atom1 = -1, min_atom2 = -1;
        for (unsigned int iatom1 = 0; iatom1 < body1.atoms().size(); iatom1++) {
            const Atom& atom1 = body1.atoms(iatom1);
            if (atom1.element != constants::symbols::carbon) {continue;}

            for (unsigned int iatom2 = 0; iatom2 < body2.atoms().size(); iatom2++) {
                const Atom& atom2 = body2.atoms(iatom2);
                if (atom2.element != constants::symbols::carbon) {continue;}

                double dist = atom1.distance(atom2);
                if (dist > min_dist) {continue;}

                min_dist = dist;
                min_atom1 = iatom1;
                min_atom2 = iatom2;
            }
        }

        rigidbody::DistanceConstraint constraint(this, ibody1, ibody2, min_atom1, min_atom2);
        add_constraint(std::move(constraint));

        if (setting::general::verbose) {
            std::cout << "\tConstraint created between bodies " << ibody1 << " and " << ibody2 << " on atoms " << body1.atoms(min_atom1).name << " and " << body2.atoms(min_atom2).name << std::endl;
        }
    }
    if (constraints.empty()) {
        throw except::unexpected("RigidBody::generate_simple_constraints: No constraints were generated. This is probably a bug.");
    }

    generate_constraint_map();
}

void RigidBody::add_constraint(std::shared_ptr<DistanceConstraint> constraint) {
    constraints.push_back(constraint);
}

void RigidBody::add_constraint(DistanceConstraint&& constraint) {
    constraints.push_back(std::make_shared<DistanceConstraint>(std::move(constraint)));
}

void RigidBody::add_constraint(unsigned int ibody1, unsigned int ibody2, unsigned int iatom1, unsigned int iatom2) {
    constraints.push_back(std::make_shared<DistanceConstraint>(this, ibody1, ibody2, iatom1, iatom2));
}

double RigidBody::chi2(fitter::HydrationFitter& fitter) const {
    std::shared_ptr<Fit> result = fitter.fit();
    return result->fval;
}

std::vector<std::shared_ptr<DistanceConstraint>> RigidBody::get_constraints() const {
    return constraints;
}

std::shared_ptr<DistanceConstraint> RigidBody::get_constraint(unsigned int index) const {
    return constraints.at(index);
}

void RigidBody::generate_constraint_map() {
    if (constraint_map.size() == bodies.size()) {return;}

    for (unsigned int i = 0; i < bodies.size(); i++) {
        constraint_map[i] = std::vector<std::shared_ptr<DistanceConstraint>>();
    }

    for (const auto& constraint : get_constraints()) {
        constraint_map.at(constraint->ibody1).push_back(constraint);
        constraint_map.at(constraint->ibody2).push_back(constraint);
    }
}

void RigidBody::apply_calibration(std::shared_ptr<fitter::Fit> calibration) {
    if (setting::general::verbose) {std::cout << "\tApplying calibration to rigid body." << std::endl;}
    this->calibration = calibration;
}

std::shared_ptr<fitter::LinearFitter> RigidBody::prepare_fitter(std::string measurement_path) {
    if (calibration == nullptr) {
        fitter::ConstrainedFitter<fitter::HydrationFitter> fitter(measurement_path, get_histogram());
        fitter.set_constraints(constraints);
        return std::make_shared<fitter::LinearFitter>(std::move(fitter));
    } else {
        auto histogram = get_histogram();
        histogram.apply_water_scaling_factor(calibration->get_parameter("c"));
        fitter::ConstrainedFitter<fitter::LinearFitter> fitter(measurement_path, std::move(histogram));
        fitter.set_constraints(constraints);
        return std::make_shared<fitter::LinearFitter>(std::move(fitter));
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
