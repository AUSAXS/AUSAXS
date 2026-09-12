// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <math/detail/Format.h>

#include <iomanip>
#include <sstream>

using namespace ausaxs;

std::string detail::format_vector(const double* data, std::size_t n) {
    std::stringstream s;
    s << "( ";
    for (std::size_t i = 0; i < n; ++i) {
        s << std::setprecision(8) << data[i] << " ";
    }
    s << ")";
    return s.str();
}

std::string detail::format_matrix(const double* data, std::size_t rows, std::size_t cols) {
    std::stringstream ss;
    for (std::size_t i = 0; i < rows; ++i) {
        ss << "\t" << std::setprecision(3);
        for (std::size_t j = 0; j < cols; ++j) {
            ss << std::setw(8) << data[i*cols + j];
        }
        ss << std::endl;
    }
    return ss.str();
}
