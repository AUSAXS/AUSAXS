#pragma once

#include <rigidbody/transform/TransformStrategy.h>

#include <array>
#include <algorithm>

namespace rigidbody {
    class RigidBody;

    /**
     * @brief RigidTransform. 
     * 
     * With this transformation strategy, everything connected to the target of the transformation will be transformed as well. 
     */
    class RigidTransform : public TransformStrategy {
        public:
            /**
             * @brief Construtor. 
             */
            RigidTransform(RigidBody* rigidbody);

            /**
             * @brief Destructor.
             */
            ~RigidTransform() override;

            /**
             * @brief Rotate all bodies connected to one side of the constraint.
             * 
             * @param rad Rotation angle in radians.
             * @param constraint The constraint. 
             */
            void rotate(double rad, Constraint& constraint);

            /**
             * @brief Translate all bodies connected to one side of the constraint. 
             * 
             * @param v The translation vector. 
             * @param constraint The constraint.
             */
            void translate(double length, Constraint& constraint) override;
    };
}