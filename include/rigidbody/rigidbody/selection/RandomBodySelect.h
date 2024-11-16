#pragma once

#include <rigidbody/selection/BodySelectStrategy.h>

#include <random>

namespace ausaxs::rigidbody {
    namespace selection {
        /**
         * @brief The next body is randomly selected, and the next constraint is randomly selected from the constraints connecting to that body.
         */
        class RandomBodySelect : public BodySelectStrategy {
            public: 
                RandomBodySelect(const RigidBody* rigidbody);
                ~RandomBodySelect() override;

                std::pair<unsigned int, int> next() override; ///< @copydoc BodySelectStrategy::next()

            private:
                std::mt19937 generator;                          // The random number generator. 
                std::uniform_int_distribution<int> distribution; // The random number distribution. 
        };
    }
}