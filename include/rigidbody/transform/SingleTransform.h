#pragma once

#include <rigidbody/transform/TransformStrategy.h>

namespace rigidbody {
    class SingleTransform : public TransformStrategy {
        public:
            SingleTransform(RigidBody* rigidbody) : TransformStrategy(rigidbody) {}
            ~SingleTransform() = default;

            void rotate(const Matrix<double>&, std::shared_ptr<Constraint> constraint) override;
            void translate(const Vector3<double>& t, std::shared_ptr<Constraint> constraint) override;
    };
}