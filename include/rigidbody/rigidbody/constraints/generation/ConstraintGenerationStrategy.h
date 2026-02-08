// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/detail/RigidbodyInternalFwd.h>
#include <utility/observer_ptr.h>

#include <vector>
#include <memory>

namespace ausaxs::rigidbody::constraints {
    class ConstraintGenerationStrategy {
        public:
            ConstraintGenerationStrategy(observer_ptr<const ConstraintManager> manager);
            virtual ~ConstraintGenerationStrategy();
            virtual std::vector<std::unique_ptr<IDistanceConstraint>> generate() const = 0;

            observer_ptr<const ConstraintManager> manager;
    };
}