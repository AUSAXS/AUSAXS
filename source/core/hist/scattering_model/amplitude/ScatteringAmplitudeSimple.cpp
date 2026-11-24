// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/scattering_model/amplitude/ScatteringAmplitudeSimple.h>

#include <algorithm>
#include <cmath>

using namespace ausaxs;
using namespace ausaxs::hist;

namespace {
    std::vector<hist::detail::CompactCoordinates> make_groups(hist::detail::CompactCoordinates&& atoms, hist::detail::CompactCoordinates&& waters) {
        std::vector<hist::detail::CompactCoordinates> groups;
        groups.reserve(2);
        groups.emplace_back(std::move(atoms));
        groups.emplace_back(std::move(waters));
        return groups;
    }
}

ScatteringAmplitudeSimple::ScatteringAmplitudeSimple(hist::detail::CompactCoordinates&& atoms, hist::detail::CompactCoordinates&& waters)
    : ScatteringAmplitude(make_groups(std::move(atoms), std::move(waters))) {}

container::Container2D<double> ScatteringAmplitudeSimple::get_excess_form_factors(const std::vector<double>& q, const ScatteringParameters& params) const {
    const int n_q = static_cast<int>(q.size());
    container::Container2D<double> f(2, n_q);
    std::ranges::transform(q, f.begin(0), [] (double qk) {return std::exp(-0.5*qk*qk);});
    std::ranges::transform(f.begin(0), f.end(0), f.begin(1), [cw = params.cw] (double f0) {return cw*f0;});
    return f;
}
