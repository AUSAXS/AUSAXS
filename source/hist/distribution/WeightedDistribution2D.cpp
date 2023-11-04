#include <hist/distribution/WeightedDistribution2D.h>
#include <hist/distribution/WeightedDistribution.h>
#include <constants/Constants.h>

#include <cmath>

using namespace hist;

void WeightedDistribution2D::add(int x, float distance, constants::axes::d_type value) {
    int i = std::round(distance*constants::axes::d_inv_width);
    index(x, i) += value;
    WeightedDistribution::entries[i].add(distance);
}