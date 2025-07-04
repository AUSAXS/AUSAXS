// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/distribution/Distribution3D.h>

using namespace ausaxs;
using namespace ausaxs::hist;

Distribution3D::Distribution3D(const WeightedDistribution3D& other) : container::Container3D<constants::axes::d_type>(other.size_x(), other.size_y(), other.size_z()) {
    for (std::size_t x = 0; x < size_x(); x++) {
        for (std::size_t y = 0; y < size_y(); y++) {
            for (std::size_t z = 0; z < size_z(); z++) {
                index(x, y, z) = other.index(x, y, z).value;
            }
        }
    }
}