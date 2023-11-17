#include <hist/distribution/WeightedDistribution1D.h>
#include <hist/distribution/WeightedDistribution.h>
#include <hist/distribution/Distribution1D.h>
#include <constants/Constants.h>

#include <cmath>

using namespace hist;

WeightedDistribution1D::WeightedDistribution1D(Distribution1D&& other) : Container1D(std::move(other)) {}

void WeightedDistribution1D::add(float distance, constants::axes::d_type value) {
    int i = std::round(distance*constants::axes::d_inv_width);
    index(i) += value;
    WeightedDistribution::entries.get()[i].add(distance);
}