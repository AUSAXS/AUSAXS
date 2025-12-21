// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/distribution/WeightedDistribution2D.h>
#include <hist/distribution/Distribution2D.h>
#include <settings/HistogramSettings.h>

using namespace ausaxs;
using namespace ausaxs::hist;

WeightedDistribution2D::WeightedDistribution2D(const Distribution2D& other) : Container2D(other.size_x(), other.size_y()) {
    for (std::size_t x = 0; x < other.size_x(); x++) {
        for (std::size_t y = 0; y < other.size_y(); y++) {
            index(x, y).value = other.index(x, y);
        }
    }
}

std::vector<double> WeightedDistribution2D::get_weighted_axis() const {
    auto d_vals = Axis(0, size_y()*settings::axes::bin_width, size_y()).as_vector(0.5);  // shift by 0.5 to get bin centers
    std::vector<double> weights(size_y());
    for (std::size_t y = 0; y < size_y(); y++) {
        unsigned int count = 0;
        for (std::size_t x = 0; x < size_x(); x++) {
            weights[y] += index(x, y).bin_center;
            count += index(x, y).count;
        }
        // If count is zero, use the default bin center; otherwise use the weighted average
        weights[y] = count ? weights[y] / count : d_vals[y];
    }
    return weights;
}