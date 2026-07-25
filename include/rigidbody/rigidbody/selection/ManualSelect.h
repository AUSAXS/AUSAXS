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
                /**
                 * @param ibody The index of the body to select on every call.
                 * @param use_constraints If false, the body is always reported as unconstrained (iconstraint == -1) so it is
                 *        transformed on its own rather than dragging its constraint group along. Required when the accompanying
                 *        parameter mask targets this body's symmetries: RigidTransform applies symmetry parameters to the
                 *        constraint group's primary body, which is not necessarily the selected one.
                 */
                ManualSelect(observer_ptr<const Rigidbody> rigidbody, unsigned int ibody, bool use_constraints = true);
                ~ManualSelect() override;

                std::pair<unsigned int, int> next() override; ///< @copydoc BodySelectStrategy::next()

            private:
                unsigned int ibody; // The index of the body to be transformed.
                bool use_constraints; // Whether the body may be transformed as part of its constraint group.
        };
    }
}