#pragma once

#include <memory>
#include <vector>

namespace rigidbody {
    class DistanceConstraint;
    class ConstraintManager;

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