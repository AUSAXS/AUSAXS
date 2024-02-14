#include <hist/distribution/WeightedDistribution2D.h>
#include <hist/distribution/Distribution2D.h>
#include <constants/Constants.h>

#include <cmath>

using namespace hist;

WeightedDistribution2D::WeightedDistribution2D(const Distribution2D& other) : Container2D(other.size_x(), other.size_y()) {
    for (std::size_t x = 0; x < other.size_x(); x++) {
        for (std::size_t y = 0; y < other.size_y(); y++) {
            index(x, y).value = other.index(x, y);
        }
    }
}

void WeightedDistribution2D::add(unsigned int x, float distance, constants::axes::d_type value) {
    int i = std::round(distance*constants::axes::d_inv_width);
    index(x, i).add(distance, value);
}

void WeightedDistribution2D::add2(unsigned int x, float distance, constants::axes::d_type value) {
    int i = std::round(distance*constants::axes::d_inv_width);
    index(x, i).add2(distance, value);
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