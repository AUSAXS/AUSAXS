/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <hist/intensity_calculator/crysol/CompositeDistanceHistogramCrysol.h>

using namespace hist;

CompositeDistanceHistogramCrysol::CompositeDistanceHistogramCrysol() = default;
CompositeDistanceHistogramCrysol::CompositeDistanceHistogramCrysol(const CompositeDistanceHistogramCrysol&) = default;
CompositeDistanceHistogramCrysol::CompositeDistanceHistogramCrysol(CompositeDistanceHistogramCrysol&&) noexcept = default;
CompositeDistanceHistogramCrysol& CompositeDistanceHistogramCrysol::operator=(CompositeDistanceHistogramCrysol&&) noexcept = default;
CompositeDistanceHistogramCrysol& CompositeDistanceHistogramCrysol::operator=(const CompositeDistanceHistogramCrysol&) = default;
CompositeDistanceHistogramCrysol::~CompositeDistanceHistogramCrysol() = default;

CompositeDistanceHistogramCrysol::CompositeDistanceHistogramCrysol(
    hist::Distribution3D&& p_aa, 
    hist::Distribution3D&& p_ax, 
    hist::Distribution3D&& p_xx, 
    hist::Distribution2D&& p_aw, 
    hist::Distribution2D&& p_wx, 
    hist::Distribution1D&& p_ww,
    hist::Distribution1D&& p_tot,
    double V
) : CompositeDistanceHistogramFFExplicitBase(std::move(p_aa), std::move(p_ax), std::move(p_xx), std::move(p_aw), std::move(p_wx), std::move(p_ww), std::move(p_tot)), average_displaced_V(V) {}

CompositeDistanceHistogramCrysol::CompositeDistanceHistogramCrysol(
    hist::Distribution3D&& p_aa, 
    hist::Distribution3D&& p_ax, 
    hist::Distribution3D&& p_xx, 
    hist::Distribution2D&& p_aw, 
    hist::Distribution2D&& p_wx, 
    hist::Distribution1D&& p_ww, 
    hist::WeightedDistribution1D&& p_tot,
    double V
) : CompositeDistanceHistogramFFExplicitBase(std::move(p_aa), std::move(p_ax), std::move(p_xx), std::move(p_aw), std::move(p_wx), std::move(p_ww), std::move(p_tot)), average_displaced_V(V) {}

void CompositeDistanceHistogramCrysol::initialize() {
    ffaa_table = form_factor::storage::atomic::get_precalculated_form_factor_table();
    ffax_table = form_factor::crysol::storage::cross::generate_table();
    ffxx_table = form_factor::crysol::storage::exv::generate_table();
}

double CompositeDistanceHistogramCrysol::exv_factor(double q) const {
    // G(q) factor from CRYSOL: https://doi.org/10.1107/S0021889895007047
    double c = constexpr_math::pow(average_displaced_V, 2./3)/(4*constants::pi);
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