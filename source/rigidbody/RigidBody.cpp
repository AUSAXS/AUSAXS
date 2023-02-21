#include <Symbols.h>
#include <rigidbody/RigidBody.h>
#include <rigidbody/Parameters.h>
#include <rigidbody/RigidTransform.h>
#include <rigidbody/SequentialSelect.h>
#include <rigidbody/SimpleParameterGeneration.h>
#include <rigidbody/RandomSelect.h>
#include <rigidbody/ConstrainedFitter.h>
#include <utility/Exceptions.h>
#include <math/Matrix.h>
#include <math/MatrixUtils.h>
#include <io/XYZWriter.h>

using namespace rigidbody;

RigidBody::RigidBody(Protein& protein) : protein(protein) {
    // Set body transformation strategy
    switch (setting::rigidbody::tsc) {
        case setting::rigidbody::RigidTransform:
            transform = std::make_unique<RigidTransform>(this); 
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
            body_selector = std::make_unique<RandomSelect>(protein);
            break;
        case setting::rigidbody::SequentialSelect:
            body_selector = std::make_unique<SequentialSelect>(protein);
            break;
        default: 
            throw except::unknown_argument("RigidBody::RigidBody: Unknown BodySelectStrategy.");
    }
}

#include <sstream>
void RigidBody::optimize(std::string measurement_path) {
    protein.generate_new_hydration();
    fitter::ConstrainedFitter fitter(measurement_path, protein.get_histogram());
    fitter.set_constraints(std::move(constraints));
    double best_chi2 = fitter.fit()->fval;

    if (setting::general::verbose) {
        utility::print_info("\nStarting rigid body optimization.");
        std::cout << "\tInitial chi2: " << best_chi2 << std::endl;
    }

    Parameters params(protein);
    std::shared_ptr<Grid> grid = protein.get_grid();
    io::XYZWriter writer(setting::general::output + "trajectory.xyz");

    auto gen_print = [] (const Protein& protein) {
        std::stringstream ss;
        for (const auto& atom : protein.atoms()) {
            ss << atom.as_pdb();
        }
        for (const auto& water: protein.waters()) {
            ss << water.as_pdb();
        }
        return ss.str();
    };
    std::string print = gen_print(protein);

    std::shared_ptr<Grid> best_grid = std::make_shared<Grid>(*grid);
    std::vector<Water> best_waters = protein.waters();
    for (int i = 0; i < 1000; i++) {
        std::stringstream iteration_out;
        iteration_out << "\nIteration " << i << std::endl;

        // select a body to be modified this iteration
        Body& body = protein.bodies.at(body_selector->next());
        Parameter param = parameter_generator->next();

        Body old_body(body);                                            // save the old body
        grid->remove(&body);                                            // remove the body from the grid
        Matrix R = matrix::rotation_matrix(param.alpha, param.beta, param.gamma);
        body.translate(param.dx);                                       // translate the body
        body.rotate(R);                                                 // rotate the body
        grid->add(&body);                                               // add the body back to the grid
        protein.generate_new_hydration(); 

        // update the body location in the fitter
        fitter.set_scattering_hist(protein.get_histogram());
        double new_chi2 = fitter.fit()->fval;

        iteration_out << "\tchi2 for new configuration: " << new_chi2 << std::endl;
        writer.write_frame(protein);

        // if the old configuration was better
        if (new_chi2 >= best_chi2) {
            body = std::move(old_body);     // restore the old body
            *grid = *best_grid;             // restore the old grid
            protein.waters() = best_waters; // restore the old waters

            // DEBUG //
            fitter.set_scattering_hist(protein.get_histogram());
            double recalc_chi2 = fitter.fit()->fval;
            iteration_out << "\trerolled changes. chi2 is now: " << recalc_chi2 << std::endl;

            protein.generate_new_hydration();
            fitter.set_scattering_hist(protein.get_histogram());
            recalc_chi2 = fitter.fit()->fval;
            iteration_out << "\tsanity check." << std::endl; 
            iteration_out << "\t\tchi2 is now: " << recalc_chi2 << std::endl;
            iteration_out << "\t\tshould be  : " << best_chi2 << std::endl;

            std::string new_print = gen_print(protein);
            if (new_print != print) {
                std::cout << iteration_out.str();
                utility::print_warning("\t\tWarning: The body has changed.");
            }

        } else {
            // accept the changes
            best_grid = std::make_shared<Grid>(*grid);
            best_waters = protein.waters();
            best_chi2 = new_chi2;
            params.update(body.uid, param);
            print = gen_print(protein);
            std::cout << "\nIteration " << i << std::endl;
            utility::print_success("\tRigidBody::optimize: Accepted changes. New best chi2: " + std::to_string(new_chi2));
        }
    }
    protein.save(setting::general::output + "optimized.pdb");
}

void RigidBody::generate_simple_constraints() {
    if (setting::general::verbose) {utility::print_info("Generating simple constraints for rigid body optimization.");}
    for (unsigned int ibody1 = 0; ibody1 < protein.bodies.size(); ibody1++) {
        for (unsigned int ibody2 = 0; ibody2 < protein.bodies.size(); ibody2++) {
            if (ibody1 == ibody2) {continue;}
            const Body& body1 = protein.body(ibody1);
            const Body& body2 = protein.body(ibody2);

            double min_dist = std::numeric_limits<double>::max();
            unsigned int min_atom1 = 0, min_atom2 = 0;
            for (unsigned int iatom1 = 0; iatom1 < body1.atoms().size(); iatom1++) {
                const Atom& atom1 = body1.atoms(iatom1);
                for (unsigned int iatom2 = 0; iatom2 < body2.atoms().size(); iatom2++) {
                    const Atom& atom2 = body2.atoms(iatom2);
                    if (atom1.name == constants::symbols::carbon && atom2.name == constants::symbols::carbon) {
                        double dist = atom1.distance(atom2);
                        if (dist < min_dist) {
                            min_dist = dist;
                            min_atom1 = iatom1;
                            min_atom2 = iatom2;
                        }
                    }
                }
            }

            // check if the bodies are close enough for a constraint to make sense
            if (min_dist > 4) {continue;} 
            rigidbody::Constraint constraint(&protein, ibody2, ibody2, min_atom1, min_atom2);
            add_constraint(std::move(constraint));

            if (setting::general::verbose) {
                std::cout << "Constraint created between bodies " << ibody1 << " and " << ibody2 << " on atoms " << body1.atoms(min_atom1).name << " and " << body2.atoms(min_atom2).name << std::endl;
            }
        }
    }
    
}

void RigidBody::add_constraint(const rigidbody::Constraint& constraint) {
    constraints.push_back(constraint);
}

void RigidBody::add_constraint(rigidbody::Constraint&& constraint) {
    constraints.push_back(std::move(constraint));
}

double RigidBody::chi2(fitter::HydrationFitter& fitter) const {
    std::shared_ptr<Fit> result = fitter.fit();
    return result->fval;
}