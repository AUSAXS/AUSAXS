// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <math/Matrix.h>

#include <vector>

namespace ausaxs::rigidbody::sequencer::detail {
    /**
     * @brief Eigen-decomposition of a real symmetric matrix.
     *
     * values[i] is the i-th eigenvalue (ascending) and vectors[i] the corresponding
     * (unit-length) eigenvector.
     */
    struct EigenResult {
        std::vector<double> values;
        std::vector<std::vector<double>> vectors;
    };

    /**
     * @brief Compute all eigenvalues and eigenvectors of a real symmetric matrix.
     *
     * Uses the cyclic Jacobi rotation method, which is robust and accurate for the small dense
     * symmetric matrices used by the symmetry fitter (up to 9x9). The input is assumed symmetric;
     * it is not modified.
     */
    EigenResult symmetric_eigen(const Matrix<double>& A);
}
