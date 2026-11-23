// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/scattering_model/ScatteringModelSimple.h>
#include <hist/scattering_model/AmplitudeHistogram.h>
#include <utility/Exceptions.h>

#include <cmath>

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
    cw = k;
    for (int i = 0; i < p.size(); ++i) {p[i] = distance_profiles.aa.index(i) + k*distance_profiles.aw.index(i) + k*k*distance_profiles.ww.index(i);}
}

namespace {
    // the exp(-q²/2) form factor of each scatterer
    void apply_gaussian_form_factor(std::vector<double>& I, const std::vector<double>& q) {
        for (std::size_t i = 0; i < I.size(); ++i) {I[i] *= std::exp(-q[i]*q[i]);}
    }
}

ScatteringProfile ScatteringModelSimple::transform_counts_gaussian(std::span<const double> counts) const {
    auto profile = transform_counts(counts, sinc_table.get_sinc_table());
    std::vector<double> I = profile.get_intensity();
    apply_gaussian_form_factor(I, profile.get_axis().as_vector());
    return {std::move(I), profile.get_axis()};
}

ScatteringProfile ScatteringModelSimple::debye_transform() const {
    return transform_counts_gaussian(p);
}

std::vector<double> ScatteringModelSimple::debye_transform_direct(const std::vector<double>& q) const {
    auto I = transform_counts(p, sinc_table_for(q), q);
    apply_gaussian_form_factor(I, q);
    return I;
}

ScatteringProfile ScatteringModelSimple::get_profile_aa() const {
    return transform_counts_gaussian(get_aa_counts().get_content());
}

ScatteringProfile ScatteringModelSimple::get_profile_aw() const {
    return transform_counts_gaussian(get_aw_counts().get_content());
}

ScatteringProfile ScatteringModelSimple::get_profile_ww() const {
    return transform_counts_gaussian(get_ww_counts().get_content());
}

std::vector<double> ScatteringModelSimple::amplitude_form_factors(bool form_factors, double q) const {
    if (form_factors) {
        throw except::invalid_argument("ScatteringModelSimple::amplitude_form_factors: The scatterers must be weighted, not resolved by form factor type.");
    }

    // the intensity is exp(-q^2) sum_ij w_i w_j sinc(q r_ij), so each scatterer is a point of weight w with the amplitude exp(-q^2/2)
    double g = std::exp(-0.5*q*q);
    std::vector<double> res(2);
    res[AmplitudeScatterers::atom_class] = g;
    res[AmplitudeScatterers::water_class] = cw*g;
    return res;
}
