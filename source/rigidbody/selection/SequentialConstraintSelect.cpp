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
    if (symmetry_only(mask)) {return next_symmetry_target(isymmetry_target);}

    unsigned int M = rigidbody->constraints->get_body_constraints(ibody).size();

    if (iconstraint == M) {
        ibody = (ibody + 1) % size_body();
        iconstraint = 0;
    }

    return {ibody, static_cast<int>(iconstraint++), -1};
}
