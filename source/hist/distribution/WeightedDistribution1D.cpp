#include <hist/distribution/WeightedDistribution1D.h>
#include <hist/distribution/WeightedDistribution.h>
#include <constants/Constants.h>

#include <cmath>

using namespace hist;

void WeightedDistribution1D::add(float distance, constants::axes::d_type value) {
    int i = std::round(distance*constants::axes::d_inv_width);
    index(i) += value;
    WeightedDistribution::entries[i].add(distance);
}