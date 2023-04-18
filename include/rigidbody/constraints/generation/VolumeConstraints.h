#pragma once

#include <rigidbody/constraints/generation/ConstraintGenerationStrategy.h>

namespace rigidbody {
    class VolumeConstraints : public ConstraintGenerationStrategy {
        public:
            /**
             * @brief Default constructor.
             */
            VolumeConstraints() = default;

            /**
             * @brief Generate a constraint.
             */
            std::shared_ptr<DistanceConstraint> generate() const;
    };
}