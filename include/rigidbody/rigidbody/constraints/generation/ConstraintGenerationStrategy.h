#pragma once

#include <rigidbody/detail/RigidbodyInternalFwd.h>

#include <vector>

namespace rigidbody::constraints {
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