// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/intensity_calculator/foxs/CompositeDistanceHistogramFoXS.h>
#include <settings/HistogramSettings.h>

using namespace ausaxs;
using namespace ausaxs::hist;

Limit CompositeDistanceHistogramFoXS::get_excluded_volume_scaling_factor_limits() const {
    return {0.8, 1.2};
}

double CompositeDistanceHistogramFoXS::exv_factor(double q, double cx) {
    constexpr double rm = 1.58;
    constexpr double c = rm*rm/(4*std::numbers::pi);
    return std::pow(cx, 3)*std::exp(-c*(std::pow(cx, 2) - 1)*q*q);
}

double CompositeDistanceHistogramFoXS::exv_factor(double q) const {
    return exv_factor(q, free_params.cx);
}