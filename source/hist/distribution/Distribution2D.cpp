#include <hist/distribution/Distribution2D.h>

#include <cmath>

using namespace hist;

Distribution2D::Distribution2D(const WeightedDistribution2D& other) : container::Container2D<constants::axes::d_type>(other.size_x(), other.size_y()) {
    for (int x = 0; x < size_x(); x++) {
        for (int y = 0; y < size_y(); y++) {
            index(x, y) = other.index(x, y).count;
        }
    }
}

void Distribution2D::add(unsigned int x, float distance, constants::axes::d_type value) {index(x, std::round(distance)) += value;}
void Distribution2D::add(unsigned int x, int32_t i, constants::axes::d_type value) {index(x, i) += value;}