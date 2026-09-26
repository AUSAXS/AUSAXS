// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/distribution/WeightedDistribution3D.h>

#include <hist/distribution/Distribution3D.h>
#include <settings/HistogramSettings.h>

#include <algorithm>
#include <cstdint>

using namespace ausaxs;
using namespace ausaxs::hist;

template<Shape S>
WeightedDistribution3D<S>::WeightedDistribution3D(const Distribution3D<S>& other) : container::Container3D<detail::WeightedEntry, S>(other.size_x(), other.size_y(), other.size_z()) {
    // both share the same layout, so the entries can be copied in storage order
    std::ranges::transform(other, *this, this->begin(), [] (double v, detail::WeightedEntry e) {e.value = v; return e;});
}

template<Shape S>
std::vector<double> WeightedDistribution3D<S>::get_weights() const {
    auto d_vals = Axis(0, this->size_z()*settings::axes::bin_width, this->size_z()).as_vector();
    std::vector<double> weights(this->size_z());
    std::vector<std::int64_t> counts(this->size_z());
    // every stored row once; for a triangular distribution that is each unordered pair once
    for (auto row : this->rows()) {
        for (int z = 0; z < this->size_z(); z++) {
            weights[z] += row[z].bin_center;
            counts[z] += row[z].count;
        }
    }

    for (int z = 0; z < this->size_z(); z++) {
        // NOLINTNEXTLINE - this is a small optimization to both avoid dividing by zero and correctly handle the case where count is zero
        weights[z] = !weights[z]*d_vals[z] + weights[z]/(!counts[z] + counts[z]);
    }
    return weights;
}

template class hist::WeightedDistribution3D<Shape::Square>;
template class hist::WeightedDistribution3D<Shape::Triangular>;
