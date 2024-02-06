#include <hist/distribution/Distribution1D.h>
#include <hist/distribution/WeightedDistribution1D.h>

#include <cmath>

using namespace hist;

Distribution1D::Distribution1D(const WeightedDistribution1D& other) : container::Container1D<constants::axes::d_type>(other.get_bins()) {}

std::vector<constants::axes::d_type> Distribution1D::as_vector() const {
    return this->data;
}

void Distribution1D::add(float distance, constants::axes::d_type value) {index(std::round(distance)) += value;}
void Distribution1D::add(int32_t i, constants::axes::d_type value) {index(i) += value;}