#include <hist/distribution/WeightedDistribution2D.h>
#include <hist/distribution/WeightedDistribution.h>
#include <hist/distribution/Distribution2D.h>
#include <constants/Constants.h>

#include <cmath>

using namespace hist;

WeightedDistribution2D::WeightedDistribution2D(Distribution2D&& other) : Container2D(std::move(other)) {}
void WeightedDistribution2D::add(int x, float distance, constants::axes::d_type value) {
    int i = std::round(distance*constants::axes::d_inv_width);
    index(x, i) += value;
    WeightedDistribution::entries.get()[i].add(distance);
}