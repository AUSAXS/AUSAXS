#include <hist/distribution/WeightedDistribution2D.h>
#include <hist/distribution/Distribution2D.h>
#include <constants/Constants.h>

#include <cmath>

using namespace hist;

WeightedDistribution2D::WeightedDistribution2D(Distribution2D& other) : Container2D(other.size_x(), other.size_y()) {
    for (int x = 0; x < other.size_x(); x++) {
        for (int y = 0; y < other.size_y(); y++) {
            index(x, y).count = other.index(x, y);
        }
    }
}

void WeightedDistribution2D::add(int x, float distance, constants::axes::d_type value) {
    int i = std::round(distance*constants::axes::d_inv_width);
    index(x, i).count += value;
    index(x, i).content += distance;
}