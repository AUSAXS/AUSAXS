// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/distribution/Distribution2D.h>

using namespace ausaxs;
using namespace ausaxs::hist;

Distribution2D::Distribution2D(const WeightedDistribution2D& other) : container::Container2D<constants::axes::d_type>(other.size_x(), other.size_y()) {
    for (std::size_t x = 0; x < size_x(); x++) {
        for (std::size_t y = 0; y < size_y(); y++) {
            index(x, y) = other.index(x, y).value;
        }
    }
}

constants::axes::d_type& Distribution2D::get_content(int i, int j) {return index(i, j);}
const constants::axes::d_type& Distribution2D::get_content(int i, int j) const {return index(i, j);}