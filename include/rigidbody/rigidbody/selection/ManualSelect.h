// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/selection/BodySelectStrategy.h>

namespace ausaxs::rigidbody {
    namespace selection {
        /**
         * @brief A fixed target is selected on every call.
         *
         * Selecting a body picks a random one of its constraints each call, as any other body strategy would. Selecting one of its declared symmetries
         * instead pins the optimizer to that single slot.
         */
        class ManualSelect : public BodySelectStrategy {
            public:
                /**
                 * @param ibody The index of the body to select on every call.
                 */
                ManualSelect(observer_ptr<const Rigidbody> rigidbody, unsigned int ibody);

                /**
                 * @brief Select one declared symmetry of a body, rather than the body itself.
                 *
                 * The target never carries a constraint: driving a symmetry leaves the host body's own atoms where they are, so there is no rigid group for
                 * the move to propagate to.
                 *
                 * @param ibody The index of the body declaring the symmetry.
                 * @param isymmetry The symmetry's slot within that body.
                 * @throws except::invalid_argument if that slot cannot be driven, i.e. it is a non-owning view onto a symmetry shared with other bodies.
                 */
                ManualSelect(observer_ptr<const Rigidbody> rigidbody, unsigned int ibody, unsigned int isymmetry);

                ~ManualSelect() override;

                Target next(const ParameterMask& mask) override; ///< @copydoc BodySelectStrategy::next()

            private:
                unsigned int ibody; // The index of the body to be transformed.
                int isymmetry = -1; // The symmetry slot to drive, or -1 to transform the body itself.
        };
    }
}
