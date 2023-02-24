#include <Symbols.h>
#include <rigidbody/RigidBody.h>
#include <rigidbody/parameters/Parameters.h>
#include <rigidbody/parameters/SimpleParameterGeneration.h>
#include <rigidbody/transform/RigidTransform.h>
#include <rigidbody/transform/SingleTransform.h>
#include <rigidbody/selection/SequentialSelect.h>
#include <rigidbody/selection/RandomSelect.h>
#include <rigidbody/selection/RandomConstraintSelect.h>
#include <rigidbody/ConstrainedFitter.h>
#include <utility/Exceptions.h>
#include <math/Matrix.h>
#include <math/MatrixUtils.h>
#include <io/XYZWriter.h>

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
        case setting::rigidbody::RigidTransform:
            transform = std::make_unique<RigidTransform>(this); 
            break;
        case setting::rigidbody::SingleTransform:
            transform = std::make_unique<SingleTransform>(this);
            break;
        default: 
            throw except::unknown_argument("RigidBody::RigidBody: Unkown TransformationStrategy.");
    }

    // Set parameter generation strategy
    switch (setting::rigidbody::pgsc) {
        case setting::rigidbody::Simple:
            parameter_generator = std::make_unique<SimpleParameterGeneration>(1000, 5, M_PI/3);
            break;
        default: 
            throw except::unknown_argument("RigidBody::RigidBody: Unknown ParameterGenerationStrategy.");
    }

    // Set body selection strategy
    switch (setting::rigidbody::bssc) {
        case setting::rigidbody::RandomSelect:
            body_selector = std::make_unique<RandomSelect>(this);
            break;
        case setting::rigidbody::RandomConstraintSelect:
            body_selector = std::make_unique<RandomConstraintSelect>(this);
            break;
        case setting::rigidbody::SequentialSelect:
            body_selector = std::make_unique<SequentialSelect>(this);
            break;
        default: 
            throw except::unknown_argument("RigidBody::RigidBody: Unknown BodySelectStrategy.");
    }
}

void RigidBody::optimize(std::string measurement_path) {
    generate_new_hydration();
    generate_constraint_map();
    auto fitter = prepare_fitter(measurement_path);
    // fitter::ConstrainedFitter<fitter::HydrationFitter> fitter(measurement_path, get_histogram());
    // fitter.set_constraints(constraints);
    double best_chi2 = fitter->fit()->fval;

    if (setting::general::verbose) {
        utility::print_info("\nStarting rigid body optimization.");
        std::cout << "\tInitial chi2: " << best_chi2 << std::endl;
    }

    Parameters params(this);
    io::XYZWriter writer(setting::general::output + "trajectory.xyz");

    // save the best configuration so we can restore it after each failed attempt
    std::shared_ptr<Grid> grid = get_grid();
    std::shared_ptr<Grid> best_grid = std::make_shared<Grid>(*grid);
    std::vector<Water> best_waters = waters();
    for (unsigned int i = 0; i < setting::rigidbody::iterations; i++) {
        writer.write_frame(this);

        // select a body to be modified this iteration
        auto [ibody, iconstraint] = body_selector->next();
        std::shared_ptr<Constraint> constraint = constraint_map.at(ibody).at(iconstraint);
        Parameter param = parameter_generator->next();

        Matrix R = matrix::rotation_matrix(param.alpha, param.beta, param.gamma);
        transform->apply(R, param.dx, constraint);
        generate_new_hydration(); 

        // update the body location in the fitter
        update_fitter(fitter);
        double new_chi2 = fitter->fit()->fval;

        writer.write_frame(this);

        // if the old configuration was better
        if (new_chi2 >= best_chi2) {
            transform->undo();          // undo the body transforms
            *grid = *best_grid;         // restore the old grid
            waters() = best_waters;     // restore the old waters

            if (i % 10 == 0 && setting::general::verbose) {
                std::cout << "\rIteration " << i << "          " << std::flush;
            }
        } else {
            // accept the changes
            best_grid = std::make_shared<Grid>(*grid);
            best_waters = waters();
            best_chi2 = new_chi2;
            // params.update(body.uid, param);
            std::cout << "\rIteration " << i << std::endl;
            utility::print_success("\tRigidBody::optimize: Accepted changes. New best chi2: " + std::to_string(new_chi2));
        }
    }

    save(setting::general::output + "optimized.pdb");
}

void RigidBody::generate_simple_constraints() {
    if (setting::general::verbose) {utility::print_info("Generating simple constraints for rigid body optimization.");}
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
            rigidbody::Constraint constraint(this, ibody1, ibody2, min_atom1, min_atom2);
            add_constraint(std::move(constraint));

            if (setting::general::verbose) {
                std::cout << "Constraint created between bodies " << ibody1 << " and " << ibody2 << " on atoms " << body1.atoms(min_atom1).name << " and " << body2.atoms(min_atom2).name << std::endl;
            }
        }
    }
    if (constraints.empty()) {
        throw except::unexpected("RigidBody::generate_simple_constraints: No constraints were generated. This is probably a bug.");
    }

    generate_constraint_map();
}

void RigidBody::add_constraint(std::shared_ptr<rigidbody::Constraint> constraint) {
    constraints.push_back(constraint);
}

void RigidBody::add_constraint(rigidbody::Constraint&& constraint) {
    constraints.push_back(std::make_shared<rigidbody::Constraint>(std::move(constraint)));
}

void RigidBody::add_constraint(unsigned int ibody1, unsigned int ibody2, unsigned int iatom1, unsigned int iatom2) {
    constraints.push_back(std::make_shared<rigidbody::Constraint>(this, ibody1, ibody2, iatom1, iatom2));
}

double RigidBody::chi2(fitter::HydrationFitter& fitter) const {
    std::shared_ptr<Fit> result = fitter.fit();
    return result->fval;
}

std::vector<std::shared_ptr<Constraint>> RigidBody::get_constraints() const {
    return constraints;
}

std::shared_ptr<Constraint> RigidBody::get_constraint(unsigned int index) const {
    return constraints.at(index);
}

void RigidBody::generate_constraint_map() {
    if (constraint_map.size() == bodies.size()) {return;}

    for (unsigned int i = 0; i < bodies.size(); i++) {
        constraint_map[i] = std::vector<std::shared_ptr<rigidbody::Constraint>>();
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
        return std::make_shared<fitter::ConstrainedFitter<fitter::HydrationFitter>>(std::move(fitter));
    } else {
        fitter::ConstrainedFitter<fitter::LinearFitter> fitter(measurement_path, get_histogram());
        fitter.set_constraints(constraints);
        return std::make_shared<fitter::ConstrainedFitter<fitter::LinearFitter>>(std::move(fitter));
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
