// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/constraints/generation/ConstraintGenerationStrategy.h>

namespace ausaxs::rigidbody::constraints {
    class LinearConstraints : public ConstraintGenerationStrategy {
        public:
            using ConstraintGenerationStrategy::ConstraintGenerationStrategy;
            std::vector<std::unique_ptr<IDistanceConstraint>> generate() const override;
    };
}