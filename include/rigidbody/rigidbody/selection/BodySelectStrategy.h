// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/RigidbodyFwd.h>
#include <rigidbody/selection/ParameterMask.h>
#include <rigidbody/selection/ParameterMaskStrategy.h>
#include <utility/observer_ptr.h>

#include <memory>
#include <utility>

namespace ausaxs::rigidbody {
    namespace selection {
        /**
         * @brief This super-class defines the interface for the body selection strategies for the rigid-body optimization. 
         * More specifically its implementations will decide in which order the bodies will be transformed by the optimization algorithm.
         */
        class BodySelectStrategy {
            public:
                BodySelectStrategy(observer_ptr<const Rigidbody> rigidbody);
                virtual ~BodySelectStrategy() = default;

                struct SelectionResult {
                    unsigned int ibody;
                    int iconstraint;
                    ParameterMask mask;
                };

                /**
                 * @brief Get the index of the next body and constraint to be transformed. 
                 * 
                 * @return A pair with the index of the body and the index of the constraint to be transformed. 
                 *         The latter is -1 if the body should be transformed independently. 
                 */
                virtual std::pair<unsigned int, int> next() = 0;

                /**
                 * @brief Like next(), but also returns a ParameterMask according to the configured mask strategy.
                 *
                 * The mask should be applied to the generated parameters before passing them to a transform strategy.
                 * The default mask strategy (AllMaskStrategy) keeps all parameters active.
                 */
                SelectionResult next_mask();

                /**
                 * @brief Replace the mask strategy used by next_masked(). Takes ownership.
                 */
                void set_mask_strategy(std::unique_ptr<ParameterMaskStrategy> strategy);

            protected: 
                observer_ptr<const Rigidbody> rigidbody;
                unsigned int N;

            private:
                std::unique_ptr<ParameterMaskStrategy> mask_strategy;
        };
    }
}