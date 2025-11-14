// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/distribution/WeightedDistribution3D.h>
#include <hist/distribution/Distribution3D.h>
#include <settings/HistogramSettings.h>

using namespace ausaxs;
using namespace ausaxs::hist;

WeightedDistribution3D::WeightedDistribution3D(const Distribution3D& other) : Container3D(other.size_x(), other.size_y(), other.size_z()) {
    // std::transform(other.begin(), other.end(), begin(), begin(), [] (const auto& val1, auto& val2) {return val2.count = val2;});
    for (std::size_t x = 0; x < other.size_x(); x++) {
        for (std::size_t y = 0; y < other.size_y(); y++) {
            for (std::size_t z = 0; z < other.size_z(); z++) {
                index(x, y, z).value = other.index(x, y, z);
            }
        }
    }
}

std::vector<double> WeightedDistribution3D::get_weights() const {
    auto d_vals = Axis(0, size_z()*settings::axes::bin_width, size_z()).as_vector();
    std::vector<double> weights(size_z());
    for (std::size_t z = 0; z < size_z(); z++) {
        unsigned int count = 0;
        for (std::size_t x = 0; x < size_x(); x++) {
            for (std::size_t y = 0; y < size_y(); y++) {
                weights[z] += index(x, y, z).bin_center;
                count += index(x, y, z).count;
            }
        }
        // this is a small optimization to both avoid dividing by zero and correctly handle the case where count is zero
        weights[z] = !weights[z]*d_vals[z] + weights[z]/(!count + count);
    }
    return weights;
}