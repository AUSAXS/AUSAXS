#include <hist/distribution/Distribution2D.h>

#include <cmath>

using namespace hist;

Distribution2D::Distribution2D(const WeightedDistribution2D& other) : container::Container2D<constants::axes::d_type>(other.size_x(), other.size_y()) {
    for (std::size_t x = 0; x < size_x(); x++) {
        for (std::size_t y = 0; y < size_y(); y++) {
            index(x, y) = other.index(x, y).value;
        }
    }
}

void Distribution2D::add(unsigned int x, float distance, constants::axes::d_type value) {index(x, std::round(distance)) += value;}
void Distribution2D::add(unsigned int x, int32_t i, constants::axes::d_type value) {index(x, i) += value;}