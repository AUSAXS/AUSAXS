// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/constraints/generation/NoConstraints.h>
#include <rigidbody/constraints/IDistanceConstraint.h>

using namespace ausaxs::rigidbody::constraints;

std::vector<std::unique_ptr<IDistanceConstraint>> NoConstraints::generate() const {return {};}