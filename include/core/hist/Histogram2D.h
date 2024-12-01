#pragma once

#include <math/Matrix.h>
#include <utility/Axis.h>

namespace ausaxs::hist {
    class Histogram2D {
        public: 
            Histogram2D() = default;

            Histogram2D(unsigned int size_x, unsigned int size_y) : data(size_x, size_y), x_axis(0, 0, size_x), y_axis(0, 0, size_y) {}

            Histogram2D(const Axis& x_axis, const Axis& y_axis) : data(x_axis.bins, y_axis.bins), x_axis(x_axis), y_axis(y_axis) {}

            std::string to_string() const;

            Matrix<double> data;
            Axis x_axis, y_axis;
    };
}