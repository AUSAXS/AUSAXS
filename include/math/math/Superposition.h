// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <math/Matrix.h>
#include <math/Vector3.h>

#include <vector>

namespace ausaxs::matrix {
    /**
     * @brief A rigid transform (apply as rotation*p + translation), with the residual RMSD it leaves.
     */
    struct SuperpositionResult {
        Matrix<double> rotation = Matrix<double>::identity(3);
        Vector3<double> translation = {0, 0, 0};
        double rmsd = 0;
    };

    /**
     * @brief The proper rotation (det = +1) best aligning two point sets with cross-covariance @p H.
     */
    Matrix<double> optimal_rotation(const Matrix<double>& H);

    /**
     * @brief The rigid transform best mapping @p from onto @p to, assuming from[i] corresponds to to[i].
     */
    SuperpositionResult superpose(const std::vector<Vector3<double>>& from, const std::vector<Vector3<double>>& to);
}
