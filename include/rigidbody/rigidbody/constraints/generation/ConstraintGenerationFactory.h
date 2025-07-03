// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <settings/RigidBodySettings.h>
#include <rigidbody/constraints/generation/ConstraintGenerationStrategy.h>

namespace ausaxs::rigidbody {
    namespace factory {
        /**
         * @brief Prepare a constraint generator. 
         */
        std::unique_ptr<constraints::ConstraintGenerationStrategy> generate_constraints(const constraints::ConstraintManager* manager);

        /**
         * @brief Prepare a constraint generator. 
         */
        std::unique_ptr<constraints::ConstraintGenerationStrategy> generate_constraints(const constraints::ConstraintManager* manager, const settings::rigidbody::ConstraintGenerationStrategyChoice& choice);
    }
}