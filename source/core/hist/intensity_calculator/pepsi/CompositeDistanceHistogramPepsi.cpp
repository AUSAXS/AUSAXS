// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/intensity_calculator/pepsi/CompositeDistanceHistogramPepsi.h>

using namespace ausaxs;

double hist::CompositeDistanceHistogramPepsi::exv_factor(double, double cx) {
    // Approximation of the G(q) factor from the Pepsi-SAXS paper, doi: 10.1107/S2059798317005745
    // This is just a Maclaurin expansion of the original expression, containing only the linear terms and no q-dependence
    // double magic_constant = 1/(4*constants::pi*constants::pi);
    double magic_constant = 1;
    double rm = 1.64;
    double c = 2*std::numbers::pi*constexpr_math::pow(4*std::numbers::pi/3, 2./3)*rm*rm*magic_constant;
    return (1 + cx*(3-c));
}

double hist::CompositeDistanceHistogramPepsi::exv_factor(double q) const {
    return exv_factor(q, free_params.cx);
}

Limit hist::CompositeDistanceHistogramPepsi::get_excluded_volume_scaling_factor_limits() const {
    return {-0.05, 0.05};
}