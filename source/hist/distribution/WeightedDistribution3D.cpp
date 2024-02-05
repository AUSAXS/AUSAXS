#include <hist/distribution/WeightedDistribution3D.h>
#include <hist/distribution/Distribution3D.h>
#include <constants/Constants.h>

#include <cmath>

using namespace hist;

WeightedDistribution3D::WeightedDistribution3D(const Distribution3D& other) : Container3D(other.size_x(), other.size_y(), other.size_z()) {
    // std::transform(other.begin(), other.end(), begin(), begin(), [] (const auto& val1, auto& val2) {return val2.count = val2;});
    for (int x = 0; x < other.size_x(); x++) {
        for (int y = 0; y < other.size_y(); y++) {
            for (int z = 0; z < other.size_z(); z++) {
                index(x, y, z).count = other.index(x, y, z);
            }
        }
    }
}

void WeightedDistribution3D::add(int x, int y, float distance, constants::axes::d_type value) {
    int i = std::round(distance*constants::axes::d_inv_width);
    index(x, y, i) += value;
    index(x, y, i).content += distance;
}