#pragma once

#include <rigidbody/constraints/generation/ConstraintGenerationStrategy.h>

namespace ausaxs::rigidbody::constraints {
    class VolumetricConstraints : public ConstraintGenerationStrategy {
        public:
            using ConstraintGenerationStrategy::ConstraintGenerationStrategy;

            /**
             * @brief Generate a constraint.
             */
            std::vector<DistanceConstraint> generate() const override;
    };
}