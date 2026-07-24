// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <math/MathFwd.h>

#include <vector>

namespace ausaxs::matrix {
    struct EigenResult {
        std::vector<double> values;                //< eigenvalues, ascending
        std::vector<std::vector<double>> vectors;  //< matching unit eigenvectors
    };

    /**
     * @brief All eigenvalues and eigenvectors of a real symmetric matrix. The input is assumed
     *        symmetric and is not modified.
     */
    EigenResult symmetric_eigen(const Matrix<double>& A);
}
