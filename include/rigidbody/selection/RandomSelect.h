#pragma once

#include <rigidbody/selection/BodySelectStrategy.h>

#include <random>

namespace rigidbody {
    /**
     * @brief Thread-safe body selection strategy. The next body is randomly selected.
     */
    class RandomSelect : public BodySelectStrategy {
        public: 
            /**
             * @brief Constructor.
             */
            RandomSelect(const RigidBody* rigidbody);

            /**
             * @brief Destructor.
             */
            ~RandomSelect() override;

            /**
             * @brief Get the index of the next body to be transformed. 
             */
            unsigned int next() override;

        private:
            std::mt19937 generator;                          // The random number generator. 
            std::uniform_int_distribution<int> distribution; // The random number distribution. 
    };
}