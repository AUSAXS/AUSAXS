// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/constraints/generation/NoConstraints.h>
#include <rigidbody/constraints/DistanceConstraint.h>

using namespace ausaxs::rigidbody::constraints;

std::vector<DistanceConstraint> NoConstraints::generate() const {return {};}