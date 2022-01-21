#include "rigidbody/RigidBody.h"

#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>

RigidBody::Parameters::Parameters(const Protein& protein) : params(protein.bodies.size()) {
    const vector<Body>& bodies = protein.bodies;
    for (unsigned int i = 0; i < params.size(); i++) {
        id_to_index[bodies[i].uid] = i;
    }
}

void RigidBody::Parameters::update(unsigned int uid, Vector3 dx, double drx, double dry, double drz) {
    params[id_to_index[uid]] = Parameter{dx, drx, dry, drz};
}

const RigidBody::Parameters::Parameter RigidBody::Parameters::get(unsigned int uid) {
    return params[id_to_index[uid]];
}

void RigidBody::optimize() {
    // ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("GSLMultiMin", "BFGS");
    // auto f = std::bind(&RigidBody::chi2, this, std::placeholders::_1);
    // ROOT::Math::Functor functor(f, 1); // declare the function to be minimized and its number of parameters
    // minimizer->SetFunction(functor);
    // minimizer->SetLimitedVariable(0, "c", 5, 1e-4, 0, 100); // scaling factor
    // minimizer->Minimize();
}

void RigidBody::driver(const string& measurement_path) {
    // Parameters params(protein);
    // double _chi2 = chi2(params);

    // IntensityFitter fitter(measurement_path, h);
    // for (unsigned int i = 0; i < 1000; i++) {
    //     // select a body to be modified this iteration
    //     Body& body = protein.bodies[body_selector->next()];


    // }
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
    Histogram h = protein.get_total_histogram();
    std::shared_ptr<Fitter::Fit> result = fitter.fit();
    return 0;
}