// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/detail/RigidbodyInternalFwd.h>

#include <vector>

namespace ausaxs::rigidbody::constraints {
    class ConstraintGenerationStrategy {
        public:
            ConstraintGenerationStrategy(const ConstraintManager* manager);
            virtual ~ConstraintGenerationStrategy();

            /**
             * @brief Generate a constraint.
             */
            virtual std::vector<DistanceConstraint> generate() const = 0;

            const ConstraintManager* manager;
    };
}