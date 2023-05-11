#pragma once

#include <rigidbody/transform/TransformStrategy.h>

namespace rigidbody {
    class SingleTransform : public TransformStrategy {
        public:
            SingleTransform(RigidBody* rigidbody);
            ~SingleTransform() override;

            /**
             * @brief Apply a transformation to the rigidbody.
             * 
             * The most recent transformation can be undone by calling undo().
             * 
             * @param M The rotation matrix.
             * @param t The translation vector. 
             * @param constraint The constraint to transform along.
             */
            void apply(const Matrix<double>& M, const Vector3<double>& t, std::shared_ptr<DistanceConstraint> constraint) override;
    };
}