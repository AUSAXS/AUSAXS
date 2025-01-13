/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <hist/distribution/WeightedDistribution2D.h>
#include <hist/distribution/Distribution2D.h>
#include <constants/Constants.h>

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
    std::vector<double> weights(size_y());
    for (std::size_t y = 0; y < size_y(); y++) {
        unsigned int count = 0;
        for (std::size_t x = 0; x < size_x(); x++) {
            weights[y] += index(x, y).bin_center;
            count += index(x, y).count;
        }
        weights[y] = !weights[y]*constants::axes::d_vals[y] + weights[y]/(!count + count); // avoid division by zero
    }
    return weights;
}