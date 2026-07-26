// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/selection/BodySelectStrategy.h>

namespace ausaxs::rigidbody {
    namespace selection {
        /**
         * @brief The next body is randomly selected, and the next constraint is randomly selected from the constraints connecting to that body.
         */
        class RandomBodySelect : public BodySelectStrategy {
            public: 
                RandomBodySelect(observer_ptr<const Rigidbody> rigidbody);
                ~RandomBodySelect() override;

                Target next(const ParameterMask& mask) override; ///< @copydoc BodySelectStrategy::next()
        };
    }
}