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
             * @brief Apply a transformation to the rigidbody.
             * 
             * The most recent transformation can be undone by calling undo().
             * 
             * @param M The rotation matrix.
             * @param t The translation vector. 
             * @param constraint The constraint to transform along.
             */
            void apply(const Matrix<double>& M, const Vector3<double>& t, std::shared_ptr<Constraint> constraint) override;
    };
}