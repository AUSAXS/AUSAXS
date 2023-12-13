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

            /**
             * @brief Apply a transformation to the rigidbody. 
             * The body to the left of the constraint (body1) is transformed, with atom2 from the right body as the pivot point.
             * 
             * The most recent transformation can be undone by calling undo().
             * 
             * @param M The rotation matrix.
             * @param t The translation vector. 
             * @param constraint The constraint to transform along.
             */
            void apply(const Matrix<double>& M, const Vector3<double>& t, constraints::DistanceConstraint& constraint) override;
    };
}