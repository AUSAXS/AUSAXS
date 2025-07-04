// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <vector>

namespace ausaxs::math {
    /**
     * @brief Interpolate points using a cubic spline. 
     *        Implementation based on some lecture notes from Aarhus University. 
     */
    class CubicSpline {
        public:
            CubicSpline(const std::vector<double>& x, const std::vector<double>& y);

            double spline(double z) const;

        private: 
            const std::vector<double> x, y;
            std::vector<double> b, c, d;

            void setup();

            /**
             * @brief Find the index @p z is supposed to be at in x.
             * @param l left bound
             * @param r right bound
             * @param z the new point on the x-axis
             */
            int search(int l, int r, double z) const;
    };
}