/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <hist/intensity_calculator/pepsi/CompositeDistanceHistogramPepsi.h>

double hist::CompositeDistanceHistogramPepsi::exv_factor(double q) const {
    // G(q) factor from CRYSOL, doi: 10.1107/S0021889895007047
    constexpr double rm = 1.62;
    constexpr double c = constexpr_math::pow(4*constants::pi/3, 3./2)*constants::pi*rm*rm;
    return std::pow(cx, 3)*std::exp(-c*(cx*cx - 1)*q*q);
}

Limit hist::CompositeDistanceHistogramPepsi::get_excluded_volume_scaling_factor_limits() const {
    return {0.865, 1.265};
}