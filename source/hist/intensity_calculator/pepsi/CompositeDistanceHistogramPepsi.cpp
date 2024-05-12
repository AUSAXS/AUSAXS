/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <hist/intensity_calculator/pepsi/CompositeDistanceHistogramPepsi.h>

double hist::CompositeDistanceHistogramPepsi::exv_factor(double q) const {
    // Approximation of the G(q) factor from the Pepsi-SAXS paper, doi: 10.1107/S2059798317005745
    // This is just a Maclaurin expansion of the original expression, containing only the linear terms and no q-dependence
    constexpr double rm = 1.64;
    constexpr double c = constexpr_math::pow(4*constants::pi/3, 2./3)*2*constants::pi*rm*rm*constants::form_factor::s_to_q_factor;
    return (1 + cx*(3-c*q*q));
}

Limit hist::CompositeDistanceHistogramPepsi::get_excluded_volume_scaling_factor_limits() const {
    return {-0.1, 0.1};
}