// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <cstddef>
#include <string>

namespace ausaxs::detail {
    /**
     * @brief Formatting helpers backing the `to_string` methods of the math containers.
     */
    std::string format_vector(const double* data, int n);

    /**
     * @brief Format a row-major matrix as one indented line per row.
     */
    std::string format_matrix(const double* data, int rows, int cols);
}