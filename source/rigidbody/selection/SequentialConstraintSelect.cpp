// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/selection/SequentialConstraintSelect.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/Rigidbody.h>

using namespace ausaxs::rigidbody::selection;

SequentialConstraintSelect::SequentialConstraintSelect(observer_ptr<const Rigidbody> rigidbody) : BodySelectStrategy(rigidbody) {}

SequentialConstraintSelect::~SequentialConstraintSelect() = default;

BodySelectStrategy::Target SequentialConstraintSelect::next(const ParameterMask& mask) {
    // a symmetry-only mask freezes the pose, so step through the drivable symmetry slots instead of the constraints; see RandomBodySelect::next
    if (symmetry_only(mask)) {
        const auto& candidates = symmetry_candidates();
        const auto& target = candidates[isymmetry_target % candidates.size()];
        isymmetry_target = (isymmetry_target + 1) % candidates.size();
        return {target.ibody, -1, static_cast<int>(target.isymmetry)};
    }

    unsigned int M = rigidbody->constraints->get_body_constraints(ibody).size();

    if (iconstraint == M) {
        ibody = (ibody + 1) % size_body();
        iconstraint = 0;
    }

    return {ibody, static_cast<int>(iconstraint++), -1};
}
