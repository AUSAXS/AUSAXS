// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/constraints/generation/ConstraintGenerationStrategy.h>

using namespace ausaxs;
using namespace rigidbody::constraints;

ConstraintGenerationStrategy::ConstraintGenerationStrategy(const ConstraintManager* manager) : manager(manager) {}
ConstraintGenerationStrategy::~ConstraintGenerationStrategy() = default;