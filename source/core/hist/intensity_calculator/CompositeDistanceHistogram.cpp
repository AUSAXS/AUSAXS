// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/intensity_calculator/CompositeDistanceHistogram.h>

using namespace ausaxs;
using namespace ausaxs::hist;

CompositeDistanceHistogram::CompositeDistanceHistogram(CompositeDistanceHistogram&&) noexcept = default;
CompositeDistanceHistogram& CompositeDistanceHistogram::operator=(CompositeDistanceHistogram&&) noexcept = default;
CompositeDistanceHistogram::~CompositeDistanceHistogram() = default;

CompositeDistanceHistogram::CompositeDistanceHistogram(
    hist::Distribution1D&& p_aa, 
    hist::Distribution1D&& p_aw, 
    hist::Distribution1D&& p_ww, 
    hist::WeightedDistribution1D&& p_tot
) : ICompositeDistanceHistogram(std::move(p_tot)), distance_profiles{.aa=std::move(p_aa), .aw=std::move(p_aw), .ww=std::move(p_ww)} {}

CompositeDistanceHistogram::CompositeDistanceHistogram(
    hist::Distribution1D&& p_aa, 
    hist::Distribution1D&& p_aw, 
    hist::Distribution1D&& p_ww, 
    hist::Distribution1D&& p_tot
) : ICompositeDistanceHistogram(std::move(p_tot)), distance_profiles{.aa=std::move(p_aa), .aw=std::move(p_aw), .ww=std::move(p_ww)} {}

const Distribution1D& CompositeDistanceHistogram::get_aa_counts() const {
    return distance_profiles.aa;
}

Distribution1D& CompositeDistanceHistogram::get_aa_counts() {
    return distance_profiles.aa;
}

const Distribution1D& CompositeDistanceHistogram::get_aw_counts() const {
    return distance_profiles.aw;
}

Distribution1D& CompositeDistanceHistogram::get_aw_counts() {
    return distance_profiles.aw;
}

const Distribution1D& CompositeDistanceHistogram::get_ww_counts() const {
    return distance_profiles.ww;
}

Distribution1D& CompositeDistanceHistogram::get_ww_counts() {
    return distance_profiles.ww;
}

void CompositeDistanceHistogram::apply_water_scaling_factor(double k) {
    for (int i = 0; i < p.size(); ++i) {p[i] = distance_profiles.aa.index(i) + k*distance_profiles.aw.index(i) + k*k*distance_profiles.ww.index(i);}
}

ScatteringProfile CompositeDistanceHistogram::get_profile_aa() const {
    return debye_sum<true>(get_aa_counts().get_content(), sinc_table.get_sinc_table());
}

ScatteringProfile CompositeDistanceHistogram::get_profile_aw() const {
    return debye_sum<true>(get_aw_counts().get_content(), sinc_table.get_sinc_table());
}

ScatteringProfile CompositeDistanceHistogram::get_profile_ww() const {
    return debye_sum<true>(get_ww_counts().get_content(), sinc_table.get_sinc_table());
}