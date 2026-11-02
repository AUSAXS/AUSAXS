// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/distribution/WeightedDistribution2D.h>

#include <hist/distribution/Distribution2D.h>

using namespace ausaxs;
using namespace ausaxs::hist;

WeightedDistribution2D::WeightedDistribution2D(const Distribution2D& other) : Container2D(other.size_x(), other.size_y()) {
    for (int x = 0; x < other.size_x(); x++) {
        for (int y = 0; y < other.size_y(); y++) {
            index(x, y).value = other.index(x, y);
        }
    }
}