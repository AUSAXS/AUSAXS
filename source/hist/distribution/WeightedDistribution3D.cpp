#include <hist/distribution/WeightedDistribution3D.h>
#include <hist/distribution/WeightedDistribution.h>
#include <constants/Constants.h>

#include <cmath>

using namespace hist;

void WeightedDistribution3D::add(int x, int y, float distance, constants::axes::d_type value) {
    int i = std::round(distance*constants::axes::d_inv_width);
    index(x, y, i) += value;
    WeightedDistribution::entries[i].add(distance);
}