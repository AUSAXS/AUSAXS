// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/intensity_calculator/crysol/CompositeDistanceHistogramCrysol.h>

using namespace ausaxs;
using namespace ausaxs::hist;

CompositeDistanceHistogramCrysol::CompositeDistanceHistogramCrysol() = default;
CompositeDistanceHistogramCrysol::CompositeDistanceHistogramCrysol(const CompositeDistanceHistogramCrysol&) = default;
CompositeDistanceHistogramCrysol::CompositeDistanceHistogramCrysol(CompositeDistanceHistogramCrysol&&) noexcept = default;
CompositeDistanceHistogramCrysol& CompositeDistanceHistogramCrysol::operator=(CompositeDistanceHistogramCrysol&&) noexcept = default;
CompositeDistanceHistogramCrysol& CompositeDistanceHistogramCrysol::operator=(const CompositeDistanceHistogramCrysol&) = default;
CompositeDistanceHistogramCrysol::~CompositeDistanceHistogramCrysol() = default;

CompositeDistanceHistogramCrysol::CompositeDistanceHistogramCrysol(
    hist::Distribution3D&& p_aa, 
    hist::Distribution2D&& p_aw, 
    hist::Distribution1D&& p_ww,
    hist::Distribution1D&& p_tot,
    double V
) : CompositeDistanceHistogramFFExplicitBase(std::move(p_aa), std::move(p_aw), std::move(p_ww), std::move(p_tot)), average_displaced_V(V) {
    initialize();
}

CompositeDistanceHistogramCrysol::CompositeDistanceHistogramCrysol(
    hist::Distribution3D&& p_aa, 
    hist::Distribution2D&& p_aw, 
    hist::Distribution1D&& p_ww, 
    hist::WeightedDistribution1D&& p_tot,
    double V
) : CompositeDistanceHistogramFFExplicitBase(std::move(p_aa), std::move(p_aw), std::move(p_ww), std::move(p_tot)), average_displaced_V(V) {
    initialize();
}

void CompositeDistanceHistogramCrysol::initialize() {
    ffaa_table = form_factor::lookup::atomic::raw::get_table();
    ffax_table = form_factor::crysol::storage::cross::generate_table();
    ffxx_table = form_factor::crysol::storage::exv::generate_table();
}

double CompositeDistanceHistogramCrysol::exv_factor(double q, double cx, double avg_displaced_V) {
    // G(q) factor from CRYSOL: https://doi.org/10.1107/S0021889895007047
    double c = constexpr_math::pow(avg_displaced_V, 2./3)/(4*std::numbers::pi);
    return std::pow(cx, 3)*std::exp(-c*(std::pow(cx, 2) - 1)*q*q);
}

double CompositeDistanceHistogramCrysol::exv_factor(double q) const {
    return exv_factor(q, free_params.cx, average_displaced_V);
}

Limit CompositeDistanceHistogramCrysol::get_excluded_volume_scaling_factor_limits() const {
    return {0.8, 1.265};
}

const form_factor::lookup::atomic::table_t& CompositeDistanceHistogramCrysol::get_ff_table() const {
    return ffaa_table;
}

const form_factor::lookup::cross::table_t& CompositeDistanceHistogramCrysol::get_ffax_table() const {
    return ffax_table;
}

const form_factor::lookup::exv::table_t& CompositeDistanceHistogramCrysol::get_ffxx_table() const {
    return ffxx_table;
}