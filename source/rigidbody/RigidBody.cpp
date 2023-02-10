#include <Symbols.h>
#include <rigidbody/RigidBody.h>
#include <rigidbody/Parameters.h>
#include <rigidbody/RigidTransform.h>
#include <rigidbody/SequentialSelect.h>
#include <rigidbody/SimpleParameterGeneration.h>
#include <rigidbody/RandomSelect.h>
#include <fitter/SimpleIntensityFitter.h>
#include <utility/Exceptions.h>
#include <math/Matrix.h>
#include <math/MatrixUtils.h>
#include <io/XYZWriter.h>

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
    SimpleIntensityFitter fitter(measurement_path, protein.get_histogram());
    double best_chi2 = fitter.fit()->fval;
    std::cout << "Initial chi2: " << best_chi2 << std::endl;

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

    for (int i = 0; i < 100; i++) {
        std::stringstream iteration_out;
        iteration_out << "\nIteration " << i << std::endl;

        // select a body to be modified this iteration
        Body& body = protein.bodies.at(body_selector->next());
        Parameter param = parameter_generator->next();

        Body old_body(body);        // save the old body
        grid->remove(&body);        // remove the body from the grid
        Matrix R = matrix::rotation_matrix(param.alpha, param.beta, param.gamma);
        body.translate(param.dx);   // translate the body
        body.rotate(R);             // rotate the body
        grid->add(&body);           // add the body back to the grid
        protein.generate_new_hydration(); 

        // update the body location in the fitter
        fitter.set_scattering_hist(protein.get_histogram());
        double new_chi2 = fitter.fit()->fval;

        iteration_out << "\tchi2 for new configuration: " << new_chi2 << std::endl;
        writer.write_frame(protein);

        // if the old configuration was better
        if (new_chi2 >= best_chi2) {
            grid->remove(&body);        // remove the body from the grid
            body = std::move(old_body); // restore the old body
            grid->add(&body);           // add the old body back to the grid
            // protein.set_grid(old_grid);
            protein.generate_new_hydration();
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
            best_chi2 = new_chi2;
            params.update(body.uid, param);
            print = gen_print(protein);
//            protein.save("temp/rigidbody/protein_" + std::to_string(i) + ".pdb");
            std::cout << "\nIteration " << i << std::endl;
            utility::print_success("\tRigidBody::optimize: Accepted changes. New best chi2: " + std::to_string(new_chi2));
        }
    }
}

void RigidBody::generate_new_hydration() {
    protein.generate_new_hydration();
}

void RigidBody::add_constraint(const Constraint& constraint) {
    constraints.push_back(constraint);
}

void RigidBody::create_constraint(const Atom* const atom1, const Atom* const atom2, const Body* const body1, const Body* const body2) {
    Constraint constraint(atom1, atom2, body1, body2);
    add_constraint(constraint);
}

void RigidBody::create_constraint(const Atom* const atom1, const Atom* const atom2) {
    auto[body1, body2] = find_host_bodies(atom1, atom2);
    create_constraint(atom1, atom2, std::move(body1), std::move(body2));
}

std::pair<const Body*, const Body*> RigidBody::find_host_bodies(const Atom* const atom1, const Atom* const atom2) const noexcept(false) {
    const Body *body1 = nullptr, *body2 = nullptr;
    const Atom a1 = *atom1; const Atom a2 = *atom2;
    for (const auto& body : protein.bodies) {
        for (const auto& atom : body.atoms()) {
            if (a1 == atom) {
                body1 = &body;
                break; // a1 and a2 *must* be from different bodies, so we break
            } else if (a2 == atom) {
                body2 = &body;
                break; // same
            }
        }
    }

    // check that both b1 and b2 were found
    if (body1 == nullptr || body2 == nullptr) {
        throw except::invalid_argument("RigidBody::create_constraint: Could not determine host bodies for the two atoms.");
    }

    return std::make_pair(body1, body2);
}

void RigidBody::create_constraint(const Atom& atom1, const Atom& atom2) {
    create_constraint(&atom1, &atom2);
}

void RigidBody::create_constraint(const Atom& atom1, const Atom& atom2, const Body& body1, const Body& body2) {
    create_constraint(&atom1, &atom2, &body1, &body2);
}

double RigidBody::chi2(IntensityFitter& fitter) const {
    std::shared_ptr<Fit> result = fitter.fit();
    return result->fval;
}