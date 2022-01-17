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

void RigidBody::create_constraint(const std::shared_ptr<Atom> const atom1, const std::shared_ptr<Atom> const atom2) {
    Constraint constraint(atom1, atom2);
    add_constraint(constraint);
}

double RigidBody::chi2() {
    // determine which body we will transform in this iteration
    Body& body = protein.bodies[body_selector->next()];

    // rotate it
    
}