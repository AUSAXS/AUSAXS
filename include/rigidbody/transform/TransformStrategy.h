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
             * @param angle The rotation angle in radians. 
             * @param constraint The constraint to rotate about. 
             */
            virtual void rotate(double rad, Constraint& constraint) = 0;

            /**
             * @brief Translate a body. 
             * 
             * @param length The distance to translate. 
             * @param constraint The constraint to translate along.
             */
            virtual void translate(double length, Constraint& constraint) = 0;

        protected: 
            RigidBody* rigidbody;
			std::unordered_map<unsigned int, std::vector<std::shared_ptr<Constraint>>> constraint_map;

            struct TransformGroup {
                std::vector<Body*> bodies;
                const Constraint* pivot;
            };

			/**
			 * @brief Generate a map of constraints for each body.
			 * 
			 * This map allows us to quickly find all constraints that apply to a given body without having to iterate over all constraints.
			 */
            void generate_constraint_map();

            /**
             * @brief Get all bodies connected by constraints to the first body of the pivot. 
             *        If we have the four bodies A - B - C - D and pivot around the BC connection, this would return the group {AB}.
             */
            TransformGroup get_connected(const Constraint& pivot);
    };
}