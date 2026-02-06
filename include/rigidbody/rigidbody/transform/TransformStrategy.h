// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/RigidbodyFwd.h>
#include <rigidbody/detail/RigidbodyInternalFwd.h>
#include <rigidbody/parameters/BodyTransformParametersRelative.h>
#include <rigidbody/parameters/BodyTransformParametersAbsolute.h>
#include <data/DataFwd.h>
#include <math/MathFwd.h>
#include <utility/observer_ptr.h>

#include <vector>

namespace ausaxs::rigidbody::transform {
    /**
     * @brief TransformStrategy. 
     * 
     * This super-class defines the interface for the body transformation strategies for the rigid-body optimization. 
     * More specifically its implementations essentially specifies how other connected bodies are affected by a transformation. 
     */
    class TransformStrategy {
        public:
            TransformStrategy(observer_ptr<Rigidbody> rigidbody);
            virtual ~TransformStrategy();

            /**
             * @brief Apply a relative transformation to a rigid group defined by a constraint.
             * 
             * The delta transformation is applied to all bodies in the smaller branch of the constraint.
             * Absolute parameters for all affected bodies are updated accordingly.
             * The most recent transformation can be undone by calling undo().
             * 
             * @param par The relative transformation to apply.
             * @param constraint The constraint to transform along.
             */
            virtual void apply(parameter::BodyTransformParametersRelative&& par, constraints::DistanceConstraint& constraint) = 0;

            /**
             * @brief Apply a relative transformation to a single unconstrained body. 
             * 
             * The delta transformation is applied to the specified body.
             * Absolute parameters for the body are updated accordingly.
             * The most recent transformation can be undone by calling undo().
             * 
             * @param par The relative transformation to apply.
             * @param ibody The index of the body to transform.
             */
            void apply(parameter::BodyTransformParametersRelative&& par, unsigned int ibody);

            /**
             * @brief Undo the previous transformation. 
             */
            virtual void undo();

        protected: 
            observer_ptr<Rigidbody> rigidbody;
            std::vector<BackupBody> bodybackup;

            /**
             * @brief Create a backup of the bodies in the group.
             */
            void backup(TransformGroup& group);

            /**
             * @brief Rotate and translate a body. 
             * 
             * @param M The rotation matrix.
             * @param t The translation vector.
             * @param group The group to apply the transformation to.
             */
            void rotate_and_translate(const Matrix<double>& M, const Vector3<double>& t, TransformGroup& group);

            /**
             * @brief Rotate a body. 
             * 
             * @param M The rotation matrix.
             * @param group The group to apply the rotation to.
             */
            virtual void rotate(const Matrix<double>& M, TransformGroup& group);

            /**
             * @brief Translate a body. 
             * 
             * @param t The translation vector. 
             * @param constraint The group to apply the translation to. 
             */
            virtual void translate(const Vector3<double>& t, TransformGroup& group);

            /**
             * @brief Set a new set of absolute symmetry parameters for a given body. 
             */
            virtual void apply_symmetry(const std::vector<symmetry::Symmetry>& symmetry, data::Body& body);
            std::vector<symmetry::Symmetry> add_symmetries(const std::vector<symmetry::Symmetry>& current, const std::vector<symmetry::Symmetry>& delta);

            /**
             * @brief Update a body with new absolute rotation and translation, handling symmetry accumulation.
             * 
             * This helper method encapsulates the common pattern:
             * 1. Get fresh body from initial_conformation
             * 2. Apply absolute rotation and translation
             * 3. Restore accumulated symmetry state
             * 4. Apply symmetry deltas (if provided)
             * 5. Update molecule and current_params
             * 
             * @param ibody The index of the body to update.
             * @param new_state The new absolute parameters for the body.
             */
            void update_body(unsigned int ibody, parameter::BodyTransformParametersAbsolute&& pars);
    };
}