// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/distribution/WeightedDistribution3D.h>

#include <hist/distribution/Distribution3D.h>

#include <algorithm>

using namespace ausaxs;
using namespace ausaxs::hist;

template<Shape S>
WeightedDistribution3D<S>::WeightedDistribution3D(const Distribution3D<S>& other) : container::Container3D<detail::WeightedEntry, S>(other.size_x(), other.size_y(), other.size_z()) {
    // both share the same layout, so the entries can be copied in storage order
    std::ranges::transform(other, *this, this->begin(), [] (double v, detail::WeightedEntry e) {e.value = v; return e;});
}

template class hist::WeightedDistribution3D<Shape::Square>;
template class hist::WeightedDistribution3D<Shape::Triangular>;
