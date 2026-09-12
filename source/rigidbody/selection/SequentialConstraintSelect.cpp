// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/selection/SequentialConstraintSelect.h>

#include <rigidbody/Rigidbody.h>
#include <rigidbody/constraints/ConstraintManager.h>

using namespace ausaxs::rigidbody::selection;

SequentialConstraintSelect::SequentialConstraintSelect(observer_ptr<const Rigidbody> rigidbody) : BodySelectStrategy(rigidbody) {}

SequentialConstraintSelect::~SequentialConstraintSelect() = default;

BodySelectStrategy::Target SequentialConstraintSelect::next(const ParameterMask& mask) {
    // a symmetry-only mask freezes the pose, so step through the drivable symmetry slots instead of the constraints; see RandomBodySelect::next
    if (symmetry_only(mask)) {return next_symmetry_target(isymmetry_target);}

    int M = static_cast<int>(rigidbody->constraints->get_body_constraints(ibody).size());

    if (iconstraint == M) {
        ibody = (ibody + 1) % size_body();
        iconstraint = 0;
    }

    return {.ibody=ibody, .iconstraint=iconstraint++, .isymmetry=-1};
}