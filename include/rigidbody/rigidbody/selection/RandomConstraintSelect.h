#pragma once

#include <rigidbody/selection/BodySelectStrategy.h>

#include <random>

namespace ausaxs::rigidbody {
    namespace selection {
        /**
         * @brief The next constraint is randomly selected, with the body being the one to which the constraint is connected.
		 *        This strategy will throw an exception if a body has no constraints.
         */
        class RandomConstraintSelect : public BodySelectStrategy {
            public: 
                RandomConstraintSelect(observer_ptr<const RigidBody> rigidbody);
                ~RandomConstraintSelect() override;

                std::pair<unsigned int, int> next() override; ///< @copydoc BodySelectStrategy::next()

            private:
                std::mt19937 generator;                          // The random number generator. 
                std::uniform_int_distribution<int> distribution; // The random number distribution. 
        };
    }
}