// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/scattering_model/ScatteringAmplitudeSimple.h>

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

std::vector<double> ScatteringAmplitudeSimple::form_factors(double q, const ScatteringParameters& params) const {
    double f = std::exp(-0.5*q*q);
    return {f, params.cw*f};
}
