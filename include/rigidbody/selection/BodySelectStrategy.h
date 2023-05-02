#pragma once

#include <rigidbody/constraints/Constraint.h>

#include <utility>

namespace rigidbody {
    class RigidBody;

    namespace selection {
        /**
         * @brief \class ConstraintSelectStrategy. 
         * 
         * This super-class defines the interface for the body selection strategies for the rigid-body optimization. 
         * More specifically its implementations will decide in which order the bodies will be transformed by the optimization algorithm.
         */
        class BodySelectStrategy {
            public:
                /**
                 * @brief Construtor. 
                 */
                BodySelectStrategy(const RigidBody* rigidbody);

                /**
                 * @brief Destructor.
                 */
                virtual ~BodySelectStrategy() = default;

                /**
                 * @brief Get the index of the next body and constraint to be transformed. 
                 */
                virtual std::pair<unsigned int, unsigned int> next() = 0;

            protected: 
                const RigidBody* rigidbody;
                unsigned int N;
        };
    }
}