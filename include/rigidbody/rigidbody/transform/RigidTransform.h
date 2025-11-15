// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/transform/TransformStrategy.h>

namespace ausaxs::rigidbody::transform {
    /**
     * @brief RigidTransform. 
     * 
     * With this transformation strategy, everything connected to the target of the transformation will be transformed as well. 
     */
    class RigidTransform : public TransformStrategy {
        public:
            RigidTransform(observer_ptr<Rigidbody> rigidbody);
            ~RigidTransform() override;

            ///< @copydoc TransformStrategy::apply(const Matrix<double>&, const Vector3<double>&, constraints::DistanceConstraint&)
            void apply(parameter::BodyTransformParameters&& par, constraints::DistanceConstraint& constraint) override;

        protected:
            /**
             * @brief Get all bodies connected by constraints to the first body of the pivot. 
             *        If we have the four bodies A - B - C - D and pivot around the BC connection, this would return the group {AB}.
             */
            TransformGroup get_connected(const constraints::DistanceConstraint& pivot);
    };
}