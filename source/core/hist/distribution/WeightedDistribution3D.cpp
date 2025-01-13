/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <hist/distribution/WeightedDistribution3D.h>
#include <hist/distribution/Distribution3D.h>
#include <constants/Constants.h>

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
    std::vector<double> weights(size_z());
    for (std::size_t z = 0; z < size_z(); z++) {
        unsigned int count = 0;
        for (std::size_t x = 0; x < size_x(); x++) {
            for (std::size_t y = 0; y < size_y(); y++) {
                weights[z] += index(x, y, z).bin_center;
                count += index(x, y, z).count;
            }
        }
        weights[z] = !weights[z]*constants::axes::d_vals[z] + weights[z]/(!count + count); // avoid division by zero
    }
    return weights;
}