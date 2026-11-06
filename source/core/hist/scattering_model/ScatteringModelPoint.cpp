// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/scattering_model/ScatteringModelPoint.h>

#include <utility>

using namespace ausaxs;
using namespace ausaxs::hist;

ScatteringModelPoint::ScatteringModelPoint(DistanceHistogram&& hist) noexcept : DistanceHistogram(std::move(hist)) {}

ScatteringProfile ScatteringModelPoint::debye_transform() const {
    return debye_sum<false>(p, sinc_table.get_sinc_table());
}

std::vector<double> ScatteringModelPoint::debye_transform_beyond_axis(const std::vector<double>& q) const {
    return debye_sum<false>(p, sinc_table_for(q), q);
}
