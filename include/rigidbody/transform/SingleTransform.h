#pragma once

#include <rigidbody/transform/TransformStrategy.h>

namespace rigidbody::transform {
    /**
     * @brief Transforms only a single body. 
     *        The body to the left of the constraint (body1) is transformed.
     */
    class SingleTransform : public TransformStrategy {
        public:
            SingleTransform(RigidBody* rigidbody);
            ~SingleTransform() override;

            ///< @copydoc TransformStrategy::apply(const Matrix<double>&, const Vector3<double>&, constraints::DistanceConstraint&)
            void apply(const Matrix<double>& M, const Vector3<double>& t, constraints::DistanceConstraint& constraint) override;

            ///< @copydoc TransformStrategy::apply(const Matrix<double>&, const Vector3<double>&, data::Body&)
            void apply(const Matrix<double>& M, const Vector3<double>& t, data::Body& body) override;
    };
}