/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include "constants/Axes.h"
#include "hist/distribution/detail/WeightedEntry.h"
#include <algorithm>
#include <hist/distribution/WeightedDistribution1D.h>
#include <hist/distribution/Distribution1D.h>
#include <constants/Constants.h>

#include <cmath>
#include <vector>
#include <cassert>

using namespace hist;

WeightedDistribution1D::WeightedDistribution1D(const Distribution1D& other) : Container1D(other.size()) {
    for (std::size_t i = 0; i < other.size(); i++) {
        index(i).value = other.index(i);
    }
}

WeightedDistribution1D::WeightedDistribution1D(const std::vector<constants::axes::d_type>& bins) : WeightedDistribution1D(Distribution1D(bins)) {}

std::vector<constants::axes::d_type> WeightedDistribution1D::as_vector() const {
    return get_content();
}

void WeightedDistribution1D::add_index(int32_t i, const detail::WeightedEntry& value) {
    index(i) += value;
}

void WeightedDistribution1D::add(float distance, constants::axes::d_type value) {
    int i = std::round(distance*constants::axes::d_inv_width);
    index(i).add(distance, value);
}

void WeightedDistribution1D::add2(float distance, constants::axes::d_type value) {
    int i = std::round(distance*constants::axes::d_inv_width);
    index(i).add2(distance, value);
}

void WeightedDistribution1D::clear(int32_t i) {
    index(i) = detail::WeightedEntry();
}

std::vector<constants::axes::d_type> WeightedDistribution1D::get_content() const {
    std::vector<constants::axes::d_type> result(size());
    for (std::size_t i = 0; i < size(); i++) {
        result[i] = index(i).value;
    }
    return result;
}

constants::axes::d_type& WeightedDistribution1D::get_content(int i) {
    return index(i).value;
}

const constants::axes::d_type& WeightedDistribution1D::get_content(int i) const {
    return index(i).value;
}

std::vector<double> WeightedDistribution1D::get_weighted_axis() const {
    Distribution1D weights(size());
    for (std::size_t i = 0; i < size(); i++) {
        weights.index(i) = (!index(i).bin_center*constants::axes::d_vals[i] + index(i).bin_center)/(!index(i).count + index(i).count); // avoid division by zero
    }
    return weights;
}

void WeightedDistribution1D::set_bin_centers(const std::vector<double>& centers) {
    assert(centers.size() == size());
    for (std::size_t i = 0; i < size(); i++) {
        index(i).bin_center = centers[i];
    }
}

WeightedDistribution1D& WeightedDistribution1D::operator+=(const WeightedDistribution1D& rhs) {
    assert(this->size() == rhs.size());
    std::transform(this->begin(), this->end(), rhs.begin(), this->begin(), std::plus<>());
    return *this;
}

WeightedDistribution1D& WeightedDistribution1D::operator-=(const WeightedDistribution1D& rhs) {
    assert(this->size() == rhs.size());
    std::transform(this->begin(), this->end(), rhs.begin(), this->begin(), std::minus<>());
    return *this;
}

WeightedDistribution1D hist::operator*(double factor, WeightedDistribution1D dist) {
    for (auto& val : dist) {
        val.value *= factor;
    }
    return dist;
}