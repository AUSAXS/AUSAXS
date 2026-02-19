// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/symmetry/Symmetry.h>
#include <math/Vector3.h>

namespace ausaxs::rigidbody::parameter {
    /**
     * @brief A small structure for storing the absolute transform parameters of a single body. 
     */
    struct BodyTransformParametersAbsolute {
        BodyTransformParametersAbsolute() : translation{0, 0, 0}, rotation{0, 0, 0} {}

        /**
         * @param dr The translation vector.
         * @param euler_angles The Euler angles.
         */
        BodyTransformParametersAbsolute(Vector3<double> dr, Vector3<double> euler_angles) : translation(std::move(dr)), rotation(std::move(euler_angles)) {}

        /**
         * @param dr The translation vector.
         * @param euler_angles The Euler angles.
         * @param symmetry_pars The symmetry parameters.
         */
        BodyTransformParametersAbsolute(Vector3<double> dr, Vector3<double> euler_angles, std::vector<symmetry::Symmetry>&& symmetry_pars)
            : translation(std::move(dr)), rotation(std::move(euler_angles)), symmetry_pars(std::move(symmetry_pars))
        {}

        void transform(const Vector3<double>& pivot, const Matrix<double>& relative_rotation) {
            translation = relative_rotation*(translation - pivot) + pivot;
            rotation = matrix::euler_angles(relative_rotation*matrix::rotation_matrix(rotation));
        }

        void transform(const Vector3<double>& relative_translation) {
            translation += relative_translation;
        }

        void transform(const Vector3<double>& pivot, const Matrix<double>& relative_rotation, const Vector3<double>& relative_translation) {
            transform(pivot, relative_rotation);
            transform(relative_translation);
        }

        bool operator==(const BodyTransformParametersAbsolute&) const = default;

        Vector3<double>                 translation;
        Vector3<double>                 rotation;
        std::vector<symmetry::Symmetry> symmetry_pars;
    };
}