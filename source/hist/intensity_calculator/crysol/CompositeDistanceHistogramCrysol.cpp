/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <hist/intensity_calculator/crysol/CompositeDistanceHistogramCrysol.h>
#include <hist/intensity_calculator/crysol/ExvFormFactorCrysol.h>
#include <settings/HistogramSettings.h>

using namespace hist;

double CompositeDistanceHistogramCrysol::exv_factor(double q) const {
    // G(q) factor from CRYSOL: https://doi.org/10.1107/S0021889895007047
    // double magic_constant = 1/(4*constants::pi*constants::pi);
    double rm = 1.62;
    double c = constexpr_math::pow(4*constants::pi/3, 3./2)*constants::pi*rm*rm;
    return std::pow(cx, 3)*std::exp(-c*(cx*cx - 1)*q*q);
}

Limit CompositeDistanceHistogramCrysol::get_excluded_volume_scaling_factor_limits() const {
    return {0.8, 1.265};
}