// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/selection/RandomConstraintSelect.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/constraints/IDistanceConstraint.h>
#include <rigidbody/Rigidbody.h>
#include <utility/Exceptions.h>
#include <utility/Random.h>

#include <random>
#include <utility>

using namespace ausaxs::rigidbody::selection;

RandomConstraintSelect::RandomConstraintSelect(observer_ptr<const Rigidbody> rigidbody) : BodySelectStrategy(rigidbody) {}

RandomConstraintSelect::~RandomConstraintSelect() = default;

BodySelectStrategy::Target RandomConstraintSelect::next(const ParameterMask& mask) {
    // a symmetry-only mask freezes the pose, so there is nothing for a constraint to propagate; draw from the drivable symmetry slots instead
    if (symmetry_only(mask)) {
        const auto& candidates = symmetry_candidates();
        std::uniform_int_distribution<std::size_t> sym_distribution(0, candidates.size()-1);
        const auto& target = candidates[sym_distribution(random::generator())];
        return {target.ibody, -1, static_cast<int>(target.isymmetry)};
    }

    // constraints connect two bodies, so deleting bodies down to a point where none remain connected empties this list
    if (rigidbody->constraints->discoverable_constraints.empty()) {
        throw except::invalid_argument("RandomConstraintSelect::next: No constraints to select from.");
    }

    // rebuild the distribution each call against the live constraint count
    std::uniform_int_distribution<int> distribution(0, rigidbody->constraints->discoverable_constraints.size()-1);
    unsigned int iconstraint = distribution(random::generator());
    const auto& constraint = rigidbody->constraints->discoverable_constraints[iconstraint];
    unsigned int ibody = constraint->ibody1;

    // find the index of the constraint in the list of constraints for the body
    for (unsigned int i = 0; i < rigidbody->constraints->get_body_constraints(ibody).size(); i++) {
        // address comparison since the DistanceConstraint comparison operator is a weak equality comparing only its contents
        if (rigidbody->constraints->get_body_constraints(ibody).at(i) == constraint.get()) {
            return {ibody, static_cast<int>(i), -1};
        }
    }
    throw except::invalid_argument("RandomConstraintSelect::next: Constraint " + std::to_string(iconstraint) + " not found");
}