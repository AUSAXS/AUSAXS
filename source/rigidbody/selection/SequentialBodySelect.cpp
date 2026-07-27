// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/selection/SequentialBodySelect.h>

using namespace ausaxs::rigidbody::selection;

SequentialBodySelect::SequentialBodySelect(observer_ptr<const Rigidbody> rigidbody) : BodySelectStrategy(rigidbody) {}

SequentialBodySelect::~SequentialBodySelect() = default;

BodySelectStrategy::Target SequentialBodySelect::next(const ParameterMask& mask) {
    // a symmetry-only mask freezes the pose, so step through the drivable symmetry slots instead of the bodies; see RandomBodySelect::next
    if (symmetry_only(mask)) {return next_symmetry_target(isymmetry_target);}

    unsigned int this_body = ibody;
    ibody = (ibody + 1) % size_body();
    return {this_body, random_constraint(this_body), -1};
}
