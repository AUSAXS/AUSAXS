#include <hist/distribution/WeightedDistribution3D.h>
#include <hist/distribution/WeightedDistribution.h>
#include <hist/distribution/Distribution3D.h>
#include <constants/Constants.h>

#include <cmath>

using namespace hist;

WeightedDistribution3D::WeightedDistribution3D(Distribution3D&& other) : Container3D(std::move(other)) {}
void WeightedDistribution3D::add(int x, int y, float distance, constants::axes::d_type value) {
    int i = std::round(distance*constants::axes::d_inv_width);
    index(x, y, i) += value;
    WeightedDistribution::entries.get()[i].add(distance);
}