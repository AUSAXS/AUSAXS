// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/constraints/generation/ConstraintGenerationStrategy.h>
#include <settings/RigidBodySettings.h>

namespace ausaxs::rigidbody::factory {
    std::unique_ptr<constraints::ConstraintGenerationStrategy> generate_constraints(
        observer_ptr<const constraints::ConstraintManager> manager, const settings::rigidbody::ConstraintGenerationStrategyChoice& choice
    );
}