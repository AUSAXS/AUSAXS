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
        BodyTransformParametersAbsolute() = default;

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

        Vector3<double>                 translation;
        Vector3<double>                 rotation;
        std::vector<symmetry::Symmetry> symmetry_pars;
    };
}