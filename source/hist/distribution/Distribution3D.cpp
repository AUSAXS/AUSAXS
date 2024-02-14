#include <hist/distribution/Distribution3D.h>

#include <cmath>

using namespace hist;

Distribution3D::Distribution3D(const WeightedDistribution3D& other) : container::Container3D<constants::axes::d_type>(other.size_x(), other.size_y(), other.size_z()) {
    for (std::size_t x = 0; x < size_x(); x++) {
        for (std::size_t y = 0; y < size_y(); y++) {
            for (std::size_t z = 0; z < size_z(); z++) {
                index(x, y, z) = other.index(x, y, z).value;
            }
        }
    }
}

void Distribution3D::add(unsigned int x, unsigned int y, float distance, constants::axes::d_type value) {index(x, y, std::round(distance)) += value;}
void Distribution3D::add2(unsigned int x, unsigned int y, float distance, constants::axes::d_type value) {index(x, y, std::round(distance)) += 2*value;}
void Distribution3D::add(unsigned int x, unsigned int y, int32_t i, constants::axes::d_type value) {index(x, y, i) += value;}