#pragma once

#include <data/Protein.h>
#include <rigidbody/Constraint.h>

namespace rigidbody {
    class RigidBody;

    /**
     * @brief TransformStrategy. 
     * 
     * This super-class defines the interface for the body transformation strategies for the rigid-body optimization. 
     * More specifically its implementations essentially specifies how other connected bodies are affected by a transformation. 
     */
    class TransformStrategy {
        public:
            /**
             * @brief Construtor. 
             */
            TransformStrategy(RigidBody* rigidbody) : rigidbody(rigidbody) {}

            /**
             * @brief Destructor.
             */
            virtual ~TransformStrategy() = default;

            /**
             * @brief Rotate a body. 
             * 
             * @param M The rotation matrix.
             */
            virtual void rotate(const Matrix<double>& M, std::shared_ptr<Constraint> constraint) = 0;

            /**
             * @brief Translate a body. 
             * 
             * @param t The translation vector. 
             * @param constraint The constraint to translate along.
             */
            virtual void translate(const Vector3<double>& t, std::shared_ptr<Constraint> constraint) = 0;

        protected: 
            RigidBody* rigidbody;

            struct TransformGroup {
                TransformGroup(std::vector<Body*> bodies, std::shared_ptr<Constraint> target, Vector3<double> pivot);
                std::vector<Body*> bodies;
                std::shared_ptr<Constraint> target;
                Vector3<double> pivot;
            };

            /**
             * @brief Get all bodies connected by constraints to the first body of the pivot. 
             *        If we have the four bodies A - B - C - D and pivot around the BC connection, this would return the group {AB}.
             */
            TransformGroup get_connected(std::shared_ptr<Constraint> pivot);
    };
}