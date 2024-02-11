#include <algorithm>
#include <hist/distribution/Distribution1D.h>
#include <hist/distribution/WeightedDistribution1D.h>

#include <cmath>

using namespace hist;

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

void Distribution1D::add(int32_t i, constants::axes::d_type value) {
    index(i) += value;
}

Distribution1D& Distribution1D::operator+=(const Distribution1D& rhs) {
    std::transform(this->begin(), this->end(), rhs.begin(), this->begin(), std::plus<>());
    return *this;
}

Distribution1D& Distribution1D::operator-=(const Distribution1D& rhs) {
    std::transform(this->begin(), this->end(), rhs.begin(), this->begin(), std::minus<>());
    return *this;
}