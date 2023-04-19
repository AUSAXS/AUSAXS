#pragma once

#include <memory>
#include <vector>

namespace rigidbody {
    class DistanceConstraint;
    class ConstraintManager;

    class ConstraintGenerationStrategy {
        public:
            ConstraintGenerationStrategy(const ConstraintManager* manager) : manager(manager) {}

            /**
             * @brief Generate a constraint.
             */
            virtual std::vector<std::shared_ptr<DistanceConstraint>> generate() const = 0;

            const ConstraintManager* manager;
    };
}