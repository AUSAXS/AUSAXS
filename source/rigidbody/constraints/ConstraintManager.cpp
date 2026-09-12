// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/constraints/ConstraintManager.h>

#include <data/Molecule.h>
#include <rigidbody/Rigidbody.h>
#include <rigidbody/constraints/DistanceConstraintBond.h>
#include <rigidbody/constraints/DistanceConstraintCM.h>
#include <rigidbody/constraints/OverlapConstraint.h>
#include <rigidbody/constraints/generation/ConstraintGenerationFactory.h>

#include <cassert>

using namespace ausaxs;
using namespace ausaxs::rigidbody::constraints;

ConstraintManager::ConstraintManager(observer_ptr<const Rigidbody> rigidbody) : molecule(&rigidbody->molecule) {
    non_discoverable_constraints.emplace_back(std::make_unique<OverlapConstraint>(molecule));
}

ConstraintManager::~ConstraintManager() = default;

const std::vector<observer_ptr<IDistanceConstraint>>& ConstraintManager::get_body_constraints(int ibody) const {
    refresh();
    assert(distance_constraints_map.contains(ibody));
    return distance_constraints_map.at(ibody);
}

void ConstraintManager::generate_constraints(settings::rigidbody::ConstraintGenerationStrategyChoice choice) {
    generate_constraints(factory::generate_constraints(this, choice));
}

void ConstraintManager::generate_constraints(std::unique_ptr<ConstraintGenerationStrategy> generator) {
    for (auto& constraint : generator->generate()) {
        discoverable_constraints.emplace_back(std::move(constraint));
    }
    invalidate();
}

void ConstraintManager::add_constraint(std::unique_ptr<Constraint> constraint) {
    if (auto* cast = dynamic_cast<DistanceConstraintCM*>(constraint.get()); cast != nullptr) {
        non_discoverable_constraints.emplace_back(std::move(constraint));
        return;
    }

    if (auto* cast = dynamic_cast<OverlapConstraint*>(constraint.get()); cast != nullptr) {
        non_discoverable_constraints.emplace_back(std::move(constraint));
        return;
    }

    auto* ptr = constraint.release();
    if (auto* cast = dynamic_cast<DistanceConstraintBond*>(ptr); cast != nullptr) {
        discoverable_constraints.emplace_back(std::unique_ptr<DistanceConstraintBond>(cast));
        invalidate();
        return;
    } if (auto* cast = dynamic_cast<DistanceConstraintAtom*>(ptr); cast != nullptr) {
        non_discoverable_constraints.emplace_back(std::unique_ptr<DistanceConstraintAtom>(cast));
        return;
    }
    delete ptr;
    throw except::invalid_argument("ConstraintManager::add_constraint: Unknown constraint type.");
}

double ConstraintManager::evaluate() const {
    double sum = 0;
    sum = std::accumulate(
        discoverable_constraints.begin(), discoverable_constraints.end(), sum, 
        [] (double sum, const std::unique_ptr<IDistanceConstraint>& constraint) {
            return sum + constraint->evaluate();
        }
    );
    sum = std::accumulate(
        non_discoverable_constraints.begin(), non_discoverable_constraints.end(), sum, 
        [] (double sum, const std::unique_ptr<Constraint>& constraint) {
            return sum + constraint->evaluate();
        }
    );
    return sum;
}

void ConstraintManager::refresh() const {
    if (!dirty) {return;}
    dirty = false;
    assert(molecule != nullptr && "ConstraintManager::refresh: Molecule must not be null.");
    distance_constraints_map.clear();
    for (int i = 0; i < molecule->size_body(); i++) {
        distance_constraints_map[i] = std::vector<observer_ptr<IDistanceConstraint>>();
    }

    for (const auto& constraint : discoverable_constraints) {
        assert(distance_constraints_map.contains(constraint->ibody1) && distance_constraints_map.contains(constraint->ibody2));
        distance_constraints_map[constraint->ibody1].emplace_back(constraint.get());
        distance_constraints_map[constraint->ibody2].emplace_back(constraint.get());
    }
}