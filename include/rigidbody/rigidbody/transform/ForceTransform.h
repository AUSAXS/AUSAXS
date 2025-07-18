// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/transform/TransformStrategy.h>

namespace ausaxs::rigidbody::transform {
    /**
     * @brief ForceTransform. 
     * 
     * This transformation strategy uses a physics-based approach to transform the rigidbody.
     * Translations are applied as forces, rotations as torques.
     */
    class ForceTransform : public TransformStrategy {
        public:
            /**
             * @brief Construtor. 
             */
            ForceTransform(RigidBody* rigidbody);

            /**
             * @brief Destructor.
             */
            ~ForceTransform() override;

            /**
             * @brief Apply a transformation to the rigidbody.
             * 
             * The most recent transformation can be undone by calling undo().
             * 
             * @param M The rotation matrix.
             * @param t The translation vector. 
             * @param constraint The constraint to transform along.
             */
            void apply(parameter::Parameter&& par, constraints::DistanceConstraint& constraint) override;
    };
}