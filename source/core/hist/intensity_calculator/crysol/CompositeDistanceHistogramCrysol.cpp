/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <hist/intensity_calculator/crysol/CompositeDistanceHistogramCrysol.h>

using namespace hist;

double CompositeDistanceHistogramCrysol::exv_factor(double q) const {
    // G(q) factor from CRYSOL: https://doi.org/10.1107/S0021889895007047
    double magic_constant = 1/(4*constants::pi*constants::pi);
    double rm = 1.62;
    double c = constexpr_math::pow(4*constants::pi/3, 3./2)*constants::pi*rm*rm*magic_constant;
    return std::pow(free_params.cx, 3)*std::exp(-c*(std::pow(free_params.cx, 2) - 1)*q*q);
}

Limit CompositeDistanceHistogramCrysol::get_excluded_volume_scaling_factor_limits() const {
    return {0.8, 1.265};
}

const form_factor::storage::atomic::table_t& CompositeDistanceHistogramCrysol::get_ff_table() const {
    return ffaa_table;
}

const form_factor::storage::cross::table_t& CompositeDistanceHistogramCrysol::get_ffax_table() const {
    return ffax_table;
}

const form_factor::storage::exv::table_t& CompositeDistanceHistogramCrysol::get_ffxx_table() const {
    return ffxx_table;
}