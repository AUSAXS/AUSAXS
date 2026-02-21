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
        BodyTransformParametersAbsolute();
        BodyTransformParametersAbsolute(const BodyTransformParametersAbsolute& other);
        BodyTransformParametersAbsolute(BodyTransformParametersAbsolute&&) noexcept;
        BodyTransformParametersAbsolute& operator=(const BodyTransformParametersAbsolute& other);
        BodyTransformParametersAbsolute& operator=(BodyTransformParametersAbsolute&&) noexcept;
        ~BodyTransformParametersAbsolute();

        /**
         * @param dr The translation vector.
         * @param euler_angles The Euler angles.
         */
        BodyTransformParametersAbsolute(Vector3<double> dr, Vector3<double> euler_angles);

        /**
         * @param dr The translation vector.
         * @param euler_angles The Euler angles.
         * @param symmetry_pars The symmetry parameters.
         */
        BodyTransformParametersAbsolute(Vector3<double> dr, Vector3<double> euler_angles, std::vector<std::unique_ptr<symmetry::ISymmetry>>&& symmetry_pars);

        void transform(const Vector3<double>& pivot, const Matrix<double>& relative_rotation);

        void transform(const Vector3<double>& relative_translation);

        void transform(const Vector3<double>& pivot, const Matrix<double>& relative_rotation, const Vector3<double>& relative_translation);

        bool operator==(const BodyTransformParametersAbsolute&) const;

        Vector3<double>                                   translation;
        Vector3<double>                                   rotation;
        std::vector<std::unique_ptr<symmetry::ISymmetry>> symmetry_pars;
    };
}