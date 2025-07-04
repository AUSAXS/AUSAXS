// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <math/Vector3.h>

namespace ausaxs::rigidbody::parameter {
    /**
     * @brief A small structure for storing a single set of parameters. 
     */
    struct Parameter {
        struct SymmetryParameter {
            Vector3<double> translation    = {0, 0, 0};
            Vector3<double> rotation_cm    = {0, 0, 0};
            Vector3<double> rotation_angle = {0, 0, 0};
        };

        Parameter() = default;

        /**
         * @brief Construct a new parameter.
         * 
         * @param dr The translation vector.
         * @param euler_angles The Euler angles.
         */
        Parameter(Vector3<double> dr, Vector3<double> euler_angles) : translation(std::move(dr)), rotation(std::move(euler_angles)) {}

        /**
         * @brief Construct a new parameter.
         * 
         * @param dr The translation vector.
         * @param euler_angles The Euler angles.
         * @param symmetry_pars The symmetry parameters.
         */
        Parameter(Vector3<double> dr, Vector3<double> euler_angles, std::vector<SymmetryParameter>&& symmetry_pars)
            : translation(std::move(dr)), rotation(std::move(euler_angles)), symmetry_pars(std::move(symmetry_pars))
        {}

        Vector3<double>                 translation;
        Vector3<double>                 rotation;
        std::vector<SymmetryParameter>  symmetry_pars;
    };
}