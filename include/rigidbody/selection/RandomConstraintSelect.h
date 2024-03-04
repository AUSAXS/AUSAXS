#pragma once

#include <rigidbody/selection/BodySelectStrategy.h>

#include <random>

namespace rigidbody {
    namespace selection {
        /**
         * @brief Thread-safe body selection strategy. The next constraint is randomly selected, with the body being the one to which the constraint is connected.
         */
        class RandomConstraintSelect : public BodySelectStrategy {
            public: 
                /**
                 * @brief Constructor.
                 */
                RandomConstraintSelect(const RigidBody* rigidbody);

                /**
                 * @brief Destructor.
                 */
                ~RandomConstraintSelect() override;

                /**
                 * @brief Get the index of the next body to be transformed. 
                 */
                std::pair<unsigned int, unsigned int> next() override;

            private:
                std::mt19937 generator;                          // The random number generator. 
                std::uniform_int_distribution<int> distribution; // The random number distribution. 
        };
    }
}