// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <math/Matrix.h>
#include <utility/Axis.h>

namespace ausaxs::hist {
    /**
     * @brief A two-dimensional histogram.
     *
     * Bin counts are stored in a matrix, together with an Axis describing the binning of each dimension.
     */
    class Histogram2D {
        public:
            Histogram2D() = default;

            /**
             * @brief Construct an empty histogram with the given number of bins along each axis.
             */
            Histogram2D(unsigned int size_x, unsigned int size_y) : data(size_x, size_y), x_axis(0, 0, size_x), y_axis(0, 0, size_y) {}

            /**
             * @brief Construct an empty histogram with the given axes.
             */
            Histogram2D(const Axis& x_axis, const Axis& y_axis) : data(x_axis.bins, y_axis.bins), x_axis(x_axis), y_axis(y_axis) {}

            std::string to_string() const;

            Matrix<double> data;
            Axis x_axis, y_axis;
    };
}