/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <hist/distribution/Distribution1D.h>
#include <hist/distribution/WeightedDistribution1D.h>

#include <cmath>
#include <algorithm>

using namespace ausaxs;
using namespace ausaxs::hist;

Distribution1D::Distribution1D(const WeightedDistribution1D& other) : container::Container1D<constants::axes::d_type>(other.get_content()) {}

std::vector<constants::axes::d_type> Distribution1D::as_vector() const {
    return this->data;
}

const std::vector<constants::axes::d_type>& Distribution1D::get_content() const {
    return this->data;
}

constants::axes::d_type& Distribution1D::get_content(int i) {
    return index(i);
}

const constants::axes::d_type& Distribution1D::get_content(int i) const {
    return index(i);
}

void Distribution1D::add(float distance, constants::axes::d_type value) {
    index(std::round(distance)) += value;
}

void Distribution1D::add2(float distance, constants::axes::d_type value) {
    index(std::round(distance)) += 2*value;
}

void Distribution1D::add_index(int32_t i, constants::axes::d_type value) {
    index(i) += value;
}

void Distribution1D::clear(int32_t i) {
    index(i) = 0;
}

Distribution1D& Distribution1D::operator+=(const Distribution1D& rhs) {
    std::transform(this->begin(), this->end(), rhs.begin(), this->begin(), std::plus<>());
    return *this;
}

Distribution1D& Distribution1D::operator-=(const Distribution1D& rhs) {
    std::transform(this->begin(), this->end(), rhs.begin(), this->begin(), std::minus<>());
    return *this;
}

Distribution1D hist::operator*(double factor, Distribution1D dist) {
    for (auto& val : dist) {
        val *= factor;
    }
    return dist;
}