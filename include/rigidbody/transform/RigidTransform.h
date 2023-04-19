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
            void apply(const Matrix<double>& M, const Vector3<double>& t, std::shared_ptr<DistanceConstraint> constraint) override;

        protected:
            /**
             * @brief Get all bodies connected by constraints to the first body of the pivot. 
             *        If we have the four bodies A - B - C - D and pivot around the BC connection, this would return the group {AB}.
             */
            TransformGroup get_connected(std::shared_ptr<DistanceConstraint> pivot);
    };
}