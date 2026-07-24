// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/selection/BodySelectStrategy.h>

namespace ausaxs::rigidbody {
    namespace selection {
        /**
         * @brief A fixed body is selected on every call, with the next constraint being a random one from the constraints connecting to that body.
         */
        class ManualSelect : public BodySelectStrategy {
            public:
                ManualSelect(observer_ptr<const Rigidbody> rigidbody, unsigned int ibody);
                ~ManualSelect() override;

                std::pair<unsigned int, int> next() override; ///< @copydoc BodySelectStrategy::next()

            private:
                unsigned int ibody; // The index of the body to be transformed.
        };
    }
}