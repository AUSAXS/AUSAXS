// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/selection/BodySelectStrategy.h>

namespace ausaxs::rigidbody {
    namespace selection {
        /**
         * @brief The next constraint is randomly selected, with the body being the one to which the constraint is connected.
		 *        This strategy will throw an exception if a body has no constraints.
         */
        class RandomConstraintSelect : public BodySelectStrategy {
            public: 
                RandomConstraintSelect(observer_ptr<const Rigidbody> rigidbody);
                ~RandomConstraintSelect() override;

                Target next(const ParameterMask& mask) override; ///< @copydoc BodySelectStrategy::next()
        };
    }
}