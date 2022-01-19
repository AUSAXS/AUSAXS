#include "rigidbody/RigidBody.h"

#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>

void RigidBody::optimize() {
    ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("GSLMultiMin", "BFGS");
    auto f = std::bind(&RigidBody::chi2, this, std::placeholders::_1);
    ROOT::Math::Functor functor(f, 1); // declare the function to be minimized and its number of parameters
    minimizer->SetFunction(functor);
    minimizer->SetLimitedVariable(0, "c", 5, 1e-4, 0, 100); // scaling factor
    minimizer->Minimize();
}

void RigidBody::add_constraint(const Constraint& constraint) {
    constraints.push_back(constraint);
}

void RigidBody::create_constraint(const std::shared_ptr<Atom> const atom1, const std::shared_ptr<Atom> const atom2, const std::shared_ptr<Body> const body1, const std::shared_ptr<Body> const body2) {
    Constraint constraint(atom1, atom2, body1, body2);
    add_constraint(constraint);
}

void RigidBody::create_constraint(const std::shared_ptr<Atom> const atom1, const std::shared_ptr<Atom> const atom2) {
    auto[body1, body2] = find_host_bodies(atom1, atom2);
    create_constraint(atom1, atom2, std::move(body1), std::move(body2));
}

std::pair<std::unique_ptr<Body>, std::unique_ptr<Body>> RigidBody::find_host_bodies(const std::shared_ptr<Atom> const atom1, const std::shared_ptr<Atom> const atom2) const noexcept(false) {
    std::unique_ptr<Body> body1 = nullptr, body2 = nullptr;
    const Atom a1 = *atom1; const Atom a2 = *atom2;
    for (const auto& body : protein.bodies) {
        for (const auto& atom : body.protein_atoms) {
            if (a1 == atom) {
                body1 = std::make_unique<Body>(body);
                break; // a1 and a2 *must* be from different bodies, so we break
            } else if (a2 == atom) {
                body2 = std::make_unique<Body>(body);
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
    std::shared_ptr<Atom> a1 = std::make_shared<Atom>(atom1);
    std::shared_ptr<Atom> a2 = std::make_shared<Atom>(atom2);
    create_constraint(a1, a2);
}

void RigidBody::create_constraint(const Atom& atom1, const Atom& atom2, const Body& body1, const Body& body2) {
    std::shared_ptr<Atom> a1 = std::make_shared<Atom>(atom1);
    std::shared_ptr<Atom> a2 = std::make_shared<Atom>(atom2);
    std::shared_ptr<Body> b1 = std::make_shared<Body>(body1);
    std::shared_ptr<Body> b2 = std::make_shared<Body>(body2);

    create_constraint(a1, a2, b1, b2);
}

double RigidBody::chi2() {
    // determine which body we will transform in this iteration
    Body& body = protein.bodies[body_selector->next()];

    // rotate it
    
}