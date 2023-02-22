#pragma once

#include <array>
#include <algorithm>

#include <data/Protein.h>
#include <rigidbody/RigidBody.h>
#include <rigidbody/transform/TransformationStrategy.h>

namespace rigidbody {
    /**
     * @brief RigidTransform. 
     * 
     * With this transformation strategy, everything connected to the target of the transformation will be transformed as well. 
     */
    class RigidTransform : public TransformationStrategy {
        public:
            /**
             * @brief Construtor. 
             */
            RigidTransform(const RigidBody* protein) : TransformationStrategy(protein) {}

            /**
             * @brief Destructor.
             */
            ~RigidTransform() override = default;

            /**
             * @brief Rotate all bodies connected to one side of the constraint.
             * 
             * @param rad Rotation angle in radians.
             * @param constraint The constraint. 
             */
            void rotate(const double rad, Constraint& constraint) override {
                std::vector<Body*> bodies = get_connected(constraint);

                Vector3 r = constraint.get_atom1().coords - constraint.get_atom2().coords;
                Vector3<double> u1, u2, u3;
                std::tie(u1, u2, u3) = r.generate_basis();

                std::for_each(bodies.begin(), bodies.end(), [&u1, &rad] (Body* body) {body->rotate(u1, rad);});
            }

            /**
             * @brief Translate all bodies connected to one side of the constraint. 
             * 
             * @param v The translation vector. 
             * @param constraint The constraint.
             */
            void translate(const double length, Constraint& constraint) override {
                std::vector<Body*> bodies = get_connected(constraint);

                Vector3 r = constraint.get_atom1().coords - constraint.get_atom2().coords;
                Vector3<double> u1, u2, u3;
                std::tie(u1, u2, u3) = r.generate_basis();

                std::for_each(bodies.begin(), bodies.end(), [&u1, &length] (Body* body) {body->translate(length*u1);});
            }
    };
}