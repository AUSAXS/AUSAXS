#include "rigidbody/RigidBody.h"
#include "rigidbody/Parameters.h"
#include "rigidbody/RigidTransform.h"
#include "rigidbody/SequentialSelect.h"
#include "rigidbody/SimpleParameterGeneration.h"
#include "rigidbody/RandomSelect.h"
#include "Exceptions.h"

#include <random>

#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>

RigidBody::RigidBody(Protein& protein) : protein(protein) {
    // Set body transformation strategy
    switch (setting::rigidbody::tsc) {
        case setting::rigidbody::RigidTransform:
            transform = std::make_unique<RigidTransform>(this); 
            break;
        default: 
            throw except::unknown_argument("Error in RigidBody::RigidBody: Unkown TransformationStrategy.");
    }

    // Set parameter generation strategy
    switch (setting::rigidbody::pgsc) {
        case setting::rigidbody::Simple:
            parameter_generator = std::make_unique<SimpleParameterGeneration>(1000, 5, M_PI/3);
            break;
        default: 
            throw except::unknown_argument("Error in RigidBody::RigidBody: Unknown ParameterGenerationStrategy.");
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
            throw except::unknown_argument("Error in RigidBody::RigidBody: Unknown BodySelectStrategy.");
    }
}

void RigidBody::optimize(const string& measurement_path) {
    generate_new_hydration();
    IntensityFitter fitter(measurement_path, protein.get_histogram());
    double _chi2 = fitter.fit()->chi2;

    std::cout << "Initial chi2: " << _chi2 << std::endl;

    Parameters params(protein);
    std::shared_ptr<Grid> grid = protein.get_grid();
    for (int i = 0; i < 100; i++) {
        std::cout << "1" << std::endl;
        // select a body to be modified this iteration
        Body& body = protein.bodies[body_selector->next()];
        Parameter param = parameter_generator->next();

        std::cout << "2" << std::endl;
        // remove the body from the grid        
        grid->remove(&body);

        // update the body to reflect the new params
        Matrix R = Matrix::rotation_matrix(param.alpha, param.beta, param.gamma);
        // body.translate(param.dx);
        // body.rotate(R);

        // add the body to the grid again
        std::cout << "3" << std::endl;
        grid->add(&body);
        protein.generate_new_hydration();

        // calculate the new chi2
        std::cout << "4" << std::endl;
        fitter.set_scattering_hist(protein.get_histogram());
        double __chi2 = fitter.fit()->chi2;

        std::cout << "chi2 for new configuration: " << __chi2 << std::endl;

        std::cout << "5" << std::endl;
        // if the old configuration was better
        if (__chi2 > _chi2) {
            std::cout << "6a" << std::endl;
            grid->remove(&body);

            Matrix R_inv = Matrix::rotation_matrix(-param.alpha, -param.beta, -param.gamma);
            // body.translate(-param.dx);
            // body.rotate(R_inv);

            grid->add(&body);
        } else {
            std::cout << "6b" << std::endl;
            // accept the changes
            _chi2 = __chi2;
            params.update(body.uid, param);
        }
        std::cout << "7" << std::endl;
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
        for (const auto& atom : body.protein_atoms) {
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
        throw except::invalid_argument("Error in RigidBody::create_constraint: Could not determine host bodies for the two atoms.");
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
    std::shared_ptr<Fitter::Fit> result = fitter.fit();
    return result->chi2;
}