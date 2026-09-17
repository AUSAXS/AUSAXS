// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/distribution/Distribution1D.h>

#include <hist/distribution/WeightedDistribution1D.h>

#include <algorithm>

using namespace ausaxs;
using namespace ausaxs::hist;

Distribution1D::Distribution1D(const WeightedDistribution1D& other) : container::Container1D<double>(other.get_content()) {}

std::vector<double> Distribution1D::as_vector() const {
    return this->data;
}

const std::vector<double>& Distribution1D::get_content() const {
    return this->data;
}

double& Distribution1D::get_content(int i) {
    return index(i);
}

const double& Distribution1D::get_content(int i) const {
    return index(i);
}

void Distribution1D::set_content(int i, double value) {
    index(i) = value;
}

void Distribution1D::clear(int32_t i) {
    index(i) = 0;
}

Distribution1D& Distribution1D::operator+=(const Distribution1D& rhs) {
    std::ranges::transform(*this, rhs, this->begin(), std::plus<>());
    return *this;
}

Distribution1D& Distribution1D::operator-=(const Distribution1D& rhs) {
    std::ranges::transform(*this, rhs, this->begin(), std::minus<>());
    return *this;
}

Distribution1D hist::operator*(double factor, Distribution1D dist) {
    for (auto& val : dist) {
        val *= factor;
    }
    return dist;
}