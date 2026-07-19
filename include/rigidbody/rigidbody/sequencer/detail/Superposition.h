// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <math/Matrix.h>
#include <math/Vector3.h>

#include <vector>

namespace ausaxs::rigidbody::sequencer::detail {
    /**
     * @brief The optimal rigid transform mapping one point set onto another.
     *
     * The transform is applied as p -> rotation*p + translation, so that
     * rotation*from[i] + translation approximates to[i]. rmsd is the residual root-mean-square
     * deviation after alignment.
     */
    struct SuperpositionResult {
        Matrix<double> rotation = Matrix<double>::identity(3);
        Vector3<double> translation = {0, 0, 0};
        double rmsd = 0;
    };

    /**
     * @brief The proper rotation R (det = +1) maximizing the Frobenius product <R, H>.
     *
     * This is the orthogonal Procrustes solution, computed via Horn's quaternion method: the
     * optimal rotation corresponds to the largest eigenvector of a symmetric 4x4 matrix built
     * from H. Because it works through quaternions it always returns a proper rotation, never a
     * reflection.
     */
    Matrix<double> optimal_rotation(const Matrix<double>& H);

    /**
     * @brief Kabsch superposition of two equally-sized, index-corresponded point sets.
     *
     * Finds the rigid transform (rotation + translation) that best maps @p from onto @p to in the
     * least-squares sense. Atom correspondence is assumed known (from[i] <-> to[i]).
     */
    SuperpositionResult superpose(const std::vector<Vector3<double>>& from, const std::vector<Vector3<double>>& to);
}
