#include <hist/distribution/WeightedDistribution1D.h>
#include <hist/distribution/Distribution1D.h>
#include <constants/Constants.h>

#include <cmath>

using namespace hist;

WeightedDistribution1D::WeightedDistribution1D(Distribution1D& other) : Container1D(other.size()) {
    for (int i = 0; i < other.size(); i++) {
        index(i).count = other.index(i);
    }
}

void WeightedDistribution1D::add(float distance, constants::axes::d_type value) {
    int i = std::round(distance*constants::axes::d_inv_width);
    index(i).count += value;
    index(i).content += distance;
}

std::vector<constants::axes::d_type> WeightedDistribution1D::get_bins() const {
    Distribution1D bins(size());
    for (int i = 0; i < size(); i++) {
        bins.index(i) = index(i).count;
    }
    return bins;
}

std::vector<double> WeightedDistribution1D::get_weights() const {
    Distribution1D weights(size());
    for (int i = 0; i < size(); i++) {
        weights.index(i) = index(i).content;
    }
    return weights;
}