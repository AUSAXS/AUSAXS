// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/RigidbodyFwd.h>
#include <rigidbody/detail/RigidbodyInternalFwd.h>
#include <rigidbody/parameters/Parameter.h>
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
            TransformStrategy(observer_ptr<RigidBody> rigidbody);
            virtual ~TransformStrategy();

            /**
             * @brief Apply a transformation to a constraint.
             * 
             * The most recent transformation can be undone by calling undo().
             * 
             * @param par The parameters to apply.
             * @param constraint The constraint to transform along.
             */
            virtual void apply(parameter::Parameter&& par, constraints::DistanceConstraint& constraint) = 0;

            /**
             * @brief Apply a transformation to a body. 
             * 
             * The most recent transformation can be undone by calling undo().
             * 
             * @param par The parameters to apply.
             * @param ibody The index of the body to transform.
             */
            virtual void apply(parameter::Parameter&& par, unsigned int ibody);

            /**
             * @brief Undo the previous transformation. 
             */
            virtual void undo();

        protected: 
            observer_ptr<RigidBody> rigidbody;
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
             * @brief Apply symmetry transformations to a body.
             */
            virtual void symmetry(std::vector<parameter::Parameter::SymmetryParameter>&& symmetry_pars, data::Body& body);
    };
}