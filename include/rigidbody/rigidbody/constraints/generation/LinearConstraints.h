#pragma once

#include <rigidbody/constraints/generation/ConstraintGenerationStrategy.h>

namespace ausaxs::rigidbody::constraints {
    class LinearConstraints : public ConstraintGenerationStrategy {
        public:
            using ConstraintGenerationStrategy::ConstraintGenerationStrategy;

            /**
             * @brief Generate a constraint.
             */
            std::vector<DistanceConstraint> generate() const override;
    };
}