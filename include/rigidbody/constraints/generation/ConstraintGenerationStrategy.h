#pragma once

#include <memory>

namespace rigidbody {
    class DistanceConstraint;
    class ConstraintManager;

    class ConstraintGenerationStrategy {
        public:
            /**
             * @brief Default constructor.
             */
            ConstraintGenerationStrategy() = default;

            /**
             * @brief Generate a constraint.
             */
            virtual std::shared_ptr<DistanceConstraint> generate() const = 0;

            ConstraintManager* manager;
    };
}