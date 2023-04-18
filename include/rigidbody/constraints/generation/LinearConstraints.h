#pragma once

#include <rigidbody/constraints/generation/ConstraintGenerationStrategy.h>

namespace rigidbody {
    class LinearConstraints : public ConstraintGenerationStrategy {
        public:
            /**
             * @brief Default constructor.
             */
            LinearConstraints() = default;

            /**
             * @brief Generate a constraint.
             */
            std::shared_ptr<DistanceConstraint> generate() const override;
    };
}