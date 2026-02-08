// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <settings/RigidBodySettings.h>
#include <rigidbody/constraints/generation/ConstraintGenerationStrategy.h>

namespace ausaxs::rigidbody::factory {
    std::unique_ptr<constraints::ConstraintGenerationStrategy> generate_constraints(observer_ptr<const constraints::ConstraintManager> manager);
    std::unique_ptr<constraints::ConstraintGenerationStrategy> generate_constraints(
        observer_ptr<const constraints::ConstraintManager> manager, const settings::rigidbody::ConstraintGenerationStrategyChoice& choice
    );
}