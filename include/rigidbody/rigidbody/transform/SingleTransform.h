// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/transform/TransformStrategy.h>

namespace ausaxs::rigidbody::transform {
    /**
     * @brief Transforms only a single body. 
     *        The body to the left of the constraint (body1) is transformed.
     */
    class SingleTransform : public TransformStrategy {
        public:
            SingleTransform(observer_ptr<RigidBody> rigidbody);
            ~SingleTransform() override;

            ///< @copydoc TransformStrategy::apply(const Matrix<double>&, const Vector3<double>&, constraints::DistanceConstraint&)
            void apply(parameter::Parameter&& par, constraints::DistanceConstraint& constraint) override;
    };
}