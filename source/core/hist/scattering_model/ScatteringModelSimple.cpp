// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/scattering_model/ScatteringModelSimple.h>

using namespace ausaxs;
using namespace ausaxs::hist;

ScatteringModelSimple::ScatteringModelSimple(ScatteringModelSimple&&) noexcept = default;
ScatteringModelSimple& ScatteringModelSimple::operator=(ScatteringModelSimple&&) noexcept = default;
ScatteringModelSimple::~ScatteringModelSimple() = default;

ScatteringModelSimple::ScatteringModelSimple(
    hist::Distribution1D&& p_aa, 
    hist::Distribution1D&& p_aw, 
    hist::Distribution1D&& p_ww, 
    hist::WeightedDistribution1D&& p_tot
) : IScatteringModel(std::move(p_tot)), distance_profiles{.aa=std::move(p_aa), .aw=std::move(p_aw), .ww=std::move(p_ww)} {}

ScatteringModelSimple::ScatteringModelSimple(
    hist::Distribution1D&& p_aa, 
    hist::Distribution1D&& p_aw, 
    hist::Distribution1D&& p_ww, 
    hist::Distribution1D&& p_tot
) : IScatteringModel(std::move(p_tot)), distance_profiles{.aa=std::move(p_aa), .aw=std::move(p_aw), .ww=std::move(p_ww)} {}

const Distribution1D& ScatteringModelSimple::get_aa_counts() const {
    return distance_profiles.aa;
}

Distribution1D& ScatteringModelSimple::get_aa_counts() {
    return distance_profiles.aa;
}

const Distribution1D& ScatteringModelSimple::get_aw_counts() const {
    return distance_profiles.aw;
}

Distribution1D& ScatteringModelSimple::get_aw_counts() {
    return distance_profiles.aw;
}

const Distribution1D& ScatteringModelSimple::get_ww_counts() const {
    return distance_profiles.ww;
}

Distribution1D& ScatteringModelSimple::get_ww_counts() {
    return distance_profiles.ww;
}

void ScatteringModelSimple::apply_water_scaling_factor(double k) {
    for (int i = 0; i < p.size(); ++i) {p[i] = distance_profiles.aa.index(i) + k*distance_profiles.aw.index(i) + k*k*distance_profiles.ww.index(i);}
}

ScatteringProfile ScatteringModelSimple::get_profile_aa() const {
    return debye_sum<true>(get_aa_counts().get_content(), sinc_table.get_sinc_table());
}

ScatteringProfile ScatteringModelSimple::get_profile_aw() const {
    return debye_sum<true>(get_aw_counts().get_content(), sinc_table.get_sinc_table());
}

ScatteringProfile ScatteringModelSimple::get_profile_ww() const {
    return debye_sum<true>(get_ww_counts().get_content(), sinc_table.get_sinc_table());
}