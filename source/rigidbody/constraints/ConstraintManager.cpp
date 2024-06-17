/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/constraints/generation/ConstraintGenerationFactory.h>
#include <rigidbody/constraints/DistanceConstraint.h>
#include <rigidbody/constraints/OverlapConstraint.h>
#include <data/record/Atom.h>
#include <data/Molecule.h>
#include <data/Body.h>

#include <cassert>

using namespace rigidbody::constraints;

ConstraintManager::ConstraintManager(data::Molecule* protein) : protein(protein), overlap_constraint(protein) {
    generate_constraints(factory::generate_constraints(this));
}

ConstraintManager::~ConstraintManager() = default;

void ConstraintManager::generate_constraints(std::unique_ptr<ConstraintGenerationStrategy> generator) {
    distance_constraints = generator->generate();
    update_constraint_map();
}

void ConstraintManager::add_constraint(std::unique_ptr<Constraint> constraint) {
    auto distance_constraint = dynamic_cast<DistanceConstraint*>(constraint.get());
    if (distance_constraint != nullptr) {
        add_constraint(std::move(*distance_constraint));
        return;
    }

    auto overlap_constraint = dynamic_cast<OverlapConstraint*>(constraint.get());
    if (overlap_constraint != nullptr) {
        add_constraint(std::move(*overlap_constraint));
        return;
    }

    throw except::invalid_argument("ConstraintManager::add_constraint: Unknown constraint type.");
}

void ConstraintManager::add_constraint(DistanceConstraint&& constraint) {
    distance_constraints.push_back(std::move(constraint));
    update_constraint_map();
}

void ConstraintManager::add_constraint(const DistanceConstraint& constraint) {
    distance_constraints.push_back(constraint);
    update_constraint_map();
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

void ConstraintManager::update_constraint_map() {
    assert(protein != nullptr && "ConstraintManager::update_constraint_map: Molecule is not set.");

    for (unsigned int i = 0; i < protein->size_body(); i++) {
        distance_constraints_map[i] = std::vector<std::reference_wrapper<DistanceConstraint>>();
    }

    for (auto& constraint : distance_constraints) {
        distance_constraints_map.at(constraint.ibody1).push_back(std::ref(constraint));
        distance_constraints_map.at(constraint.ibody2).push_back(std::ref(constraint));
    }
}