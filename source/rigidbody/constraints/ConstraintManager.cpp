#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/constraints/generation/ConstraintGenerationFactory.h>
#include <rigidbody/constraints/DistanceConstraint.h>
#include <rigidbody/constraints/OverlapConstraint.h>
#include <data/Protein.h>
#include <data/Body.h>
#include <data/Atom.h>

using namespace rigidbody;

ConstraintManager::ConstraintManager() = default;

ConstraintManager::ConstraintManager(Protein* protein) : protein(protein), overlap_constraint(protein) {
    auto generator = factory::generate_constraints(this);
    distance_constraints = generator->generate();
    generate_constraint_map();
}

ConstraintManager::~ConstraintManager() = default;

double ConstraintManager::evaluate() const {
    double chi2 = 0.0;
    for (const auto& constraint : distance_constraints) {
        chi2 += constraint.evaluate();
    }
    chi2 += overlap_constraint.evaluate();
    return chi2;
}

void ConstraintManager::generate_constraint_map() {
    if (distance_constraints_map.size() == protein->bodies.size()) {return;}

    for (unsigned int i = 0; i < protein->bodies.size(); i++) {
        distance_constraints_map[i] = std::vector<DistanceConstraint*>();
    }

    for (auto& constraint : distance_constraints) {
        distance_constraints_map.at(constraint.ibody1).push_back(&constraint);
        distance_constraints_map.at(constraint.ibody2).push_back(&constraint);
    }
}