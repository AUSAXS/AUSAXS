// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/RigidbodyFwd.h>
#include <utility/observer_ptr.h>

#include <utility>

namespace ausaxs::rigidbody {
    namespace selection {
        /**
         * @brief This super-class defines the interface for the body selection strategies for the rigid-body optimization. 
         * More specifically its implementations will decide in which order the bodies will be transformed by the optimization algorithm.
         */
        class BodySelectStrategy {
            public:
                BodySelectStrategy(observer_ptr<const RigidBody> rigidbody);
                virtual ~BodySelectStrategy() = default;

                /**
                 * @brief Get the index of the next body and constraint to be transformed. 
                 * 
                 * @return A pair with the index of the body and the index of the constraint to be transformed. 
                 *         The latter is -1 if the body should be transformed independently. 
                 */
                virtual std::pair<unsigned int, int> next() = 0;

            protected: 
                const RigidBody* rigidbody;
                unsigned int N;
        };
    }
}