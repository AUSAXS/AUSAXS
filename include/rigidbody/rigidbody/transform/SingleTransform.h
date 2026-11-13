// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/transform/TransformStrategy.h>

namespace ausaxs::rigidbody::transform {
    /**
     * @brief Transforms only a single body. 
     *        The selected body is transformed, pivoting about the constraint's atom on the other body; the rest of the system is left alone.
     */
    class SingleTransform : public TransformStrategy {
        public:
            SingleTransform(observer_ptr<Rigidbody> rigidbody);
            ~SingleTransform() override;

            ///< @copydoc TransformStrategy::apply
            void apply(const parameter::BodyTransformParametersRelative& par, observer_ptr<const constraints::IDistanceConstraint> constraint, int ibody) override;
    };
}