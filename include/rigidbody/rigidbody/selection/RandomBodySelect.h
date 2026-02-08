// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

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
                RandomBodySelect(observer_ptr<const Rigidbody> rigidbody);
                ~RandomBodySelect() override;

                std::pair<unsigned int, int> next() override; ///< @copydoc BodySelectStrategy::next()

            private:
                std::uniform_int_distribution<int> distribution; // The random number distribution. 
        };
    }
}