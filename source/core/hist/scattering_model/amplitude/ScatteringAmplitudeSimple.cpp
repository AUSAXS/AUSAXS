// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/scattering_model/amplitude/ScatteringAmplitudeSimple.h>

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
    for (int k = 0; k < n_q; ++k) {f(0, k) = std::exp(-0.5*q[k]*q[k]);}
    for (int k = 0; k < n_q; ++k) {f(1, k) = params.cw*f(0, k);}
    return f;
}
