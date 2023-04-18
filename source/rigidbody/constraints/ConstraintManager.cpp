#include <rigidbody/constraints/ConstraintManager.h>

using namespace rigidbody;

ConstraintManager::ConstraintManager(Protein* protein) : protein(protein) {}

void ConstraintManager::add_constraint(std::shared_ptr<DistanceConstraint> constraint) {
    distance_constraints.push_back(std::move(constraint));
}

void ConstraintManager::add_constraint(std::shared_ptr<OverlapConstraint> constraint) {
    overlap_constraint = std::move(constraint);
}

double ConstraintManager::evaluate() const {
    double chi2 = 0.0;
    for (const auto& constraint : distance_constraints) {
        chi2 += constraint->evaluate();
    }
    chi2 += overlap_constraint->evaluate();
    return chi2;
}