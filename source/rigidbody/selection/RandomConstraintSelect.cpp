// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/selection/RandomConstraintSelect.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/constraints/DistanceConstraint.h>
#include <rigidbody/Rigidbody.h>
#include <utility/Exceptions.h>
#include <utility/Random.h>

#include <utility>

using namespace ausaxs::rigidbody::selection;

RandomConstraintSelect::RandomConstraintSelect(observer_ptr<const Rigidbody> rigidbody) : BodySelectStrategy(rigidbody) {
    unsigned int M = rigidbody->constraints->discoverable_constraints.size();
    distribution = std::uniform_int_distribution<int>(0, M-1);
}

RandomConstraintSelect::~RandomConstraintSelect() = default;

std::pair<unsigned int, int> RandomConstraintSelect::next() {
    unsigned int iconstraint = distribution(random::generator());
    const auto& constraint = rigidbody->constraints->discoverable_constraints[iconstraint];
    unsigned int ibody = constraint->ibody1;

    // find the index of the constraint in the list of constraints for the body
    for (unsigned int i = 0; i < rigidbody->constraints->get_body_constraints(ibody).size(); i++) {
        // address comparison since the DistanceConstraint comparison operator is a weak equality comparing only its contents
        if (rigidbody->constraints->get_body_constraints(ibody).at(i) == constraint.get()) {
            return std::make_pair(ibody, i);
        }
    }
    throw except::invalid_argument("RandomConstraintSelect::next: Constraint " + std::to_string(iconstraint) + " not found");
}