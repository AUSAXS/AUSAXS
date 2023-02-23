#pragma once

#include <rigidbody/transform/TransformStrategy.h>

namespace rigidbody {
    class SingleTransform : public TransformStrategy {
        public:
            SingleTransform(RigidBody* rigidbody) : TransformStrategy(rigidbody) {}
            ~SingleTransform() = default;

            void rotate(double rad, Constraint& constraint) override;
            void translate(double length, Constraint& constraint) override;
    };
}