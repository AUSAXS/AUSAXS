// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/distribution/Distribution3D.h>

#include <algorithm>

using namespace ausaxs;
using namespace ausaxs::hist;

template<Shape S>
Distribution3D<S>::Distribution3D(const WeightedDistribution3D<S>& other) : container::Container3D<double, S>(other.size_x(), other.size_y(), other.size_z()) {
    // both share the same layout, so the entries can be copied in storage order
    std::ranges::transform(other, this->begin(), [] (const detail::WeightedEntry& e) {return e.value;});
}

template class hist::Distribution3D<Shape::Square>;
template class hist::Distribution3D<Shape::Triangular>;
