/*
This software is distributed under the GNU General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/constraints/generation/ConstraintGenerationFactory.h>
#include <rigidbody/constraints/DistanceConstraint.h>
#include <rigidbody/constraints/OverlapConstraint.h>
#include <data/record/Atom.h>
#include <data/Molecule.h>
#include <data/Body.h>

using namespace rigidbody::constraints;

ConstraintManager::ConstraintManager(data::Molecule* protein) : protein(protein), overlap_constraint(protein) {
    auto generator = factory::generate_constraints(this);
    distance_constraints = generator->generate();
    generate_constraint_map();
}

ConstraintManager::~ConstraintManager() = default;


void ConstraintManager::add_constraint(DistanceConstraint&& constraint) {
    distance_constraints.push_back(std::move(constraint));
    generate_constraint_map();
}

void ConstraintManager::add_constraint(const DistanceConstraint& constraint) {
    distance_constraints.push_back(constraint);
    generate_constraint_map();
}

void ConstraintManager::add_constraint(OverlapConstraint&& constraint) {
    overlap_constraint = std::move(constraint);
}

void ConstraintManager::add_constraint(const OverlapConstraint& constraint) {
    overlap_constraint = constraint;
}

double ConstraintManager::evaluate() const {
    double chi2 = 0.0;
    for (const auto& constraint : distance_constraints) {
        chi2 += constraint.evaluate();
    }
    chi2 += overlap_constraint.evaluate();
    return chi2;
}

void ConstraintManager::generate_constraint_map() {
    #ifdef DEBUG
        if (protein == nullptr) [[unlikely]] {throw except::nullptr_error("ConstraintManager::generate_constraint_map: Protein is not set.");}
    #endif

    for (unsigned int i = 0; i < protein->body_size(); i++) {
        distance_constraints_map[i] = std::vector<std::reference_wrapper<DistanceConstraint>>();
    }

    for (auto& constraint : distance_constraints) {
        distance_constraints_map.at(constraint.ibody1).push_back(std::ref(constraint));
        distance_constraints_map.at(constraint.ibody2).push_back(std::ref(constraint));
    }
}