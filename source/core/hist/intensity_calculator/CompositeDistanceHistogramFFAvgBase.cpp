/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <hist/intensity_calculator/CompositeDistanceHistogramFFAvgBase.h>
#include <hist/Histogram.h>
#include <table/ArrayDebyeTable.h>
#include <form_factor/FormFactor.h>
#include <form_factor/PrecalculatedFormFactorProduct.h>
#include <settings/HistogramSettings.h>

using namespace hist;

template<typename FormFactorTableType>
CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::CompositeDistanceHistogramFFAvgBase() = default;

template<typename FormFactorTableType>
CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::CompositeDistanceHistogramFFAvgBase(
    hist::Distribution3D&& p_aa, 
    hist::Distribution2D&& p_aw, 
    hist::Distribution1D&& p_ww,
    hist::Distribution1D&& p_tot
) : ICompositeDistanceHistogramExv(std::move(p_tot)), distance_profiles{.aa=std::move(p_aa), .aw=std::move(p_aw), .ww=std::move(p_ww)} {}

template<typename FormFactorTableType>
CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::CompositeDistanceHistogramFFAvgBase(
    hist::Distribution3D&& p_aa, 
    hist::Distribution2D&& p_aw, 
    hist::Distribution1D&& p_ww, 
    hist::WeightedDistribution1D&& p_tot
) : ICompositeDistanceHistogramExv(std::move(p_tot)), distance_profiles{.aa=std::move(p_aa), .aw=std::move(p_aw), .ww=std::move(p_ww)} {}

template<typename FormFactorTableType>
CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::~CompositeDistanceHistogramFFAvgBase() = default;

template<typename FormFactorTableType>
double CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::exv_factor(double) const {
    return free_params.cx;
}

template<typename FormFactorTableType>
ScatteringProfile CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::debye_transform() const {
    const auto& ff_table = get_ff_table();
    auto sinqd_table = get_sinc_table();

    // calculate the Debye scattering intensity
    Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
    unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

    std::vector<double> Iq(debye_axis.bins, 0);
    for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
        double cx = exv_factor(q);
        for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
            // atom-atom
            for (unsigned int ff2 = 0; ff2 < form_factor::get_count_without_excluded_volume(); ++ff2) {
                double aa_sum = std::inner_product(distance_profiles.aa.begin(ff1, ff2), distance_profiles.aa.end(ff1, ff2), sinqd_table->begin(q), 0.0);
                Iq[q-q0] += aa_sum*ff_table.index(ff1, ff2).evaluate(q);
            }

            // atom-exv
            double ax_sum = std::inner_product(distance_profiles.aa.begin(ff1, form_factor::exv_bin), distance_profiles.aa.end(ff1, form_factor::exv_bin), sinqd_table->begin(q), 0.0);
            Iq[q-q0] -= 2*cx*ax_sum*ff_table.index(ff1, form_factor::exv_bin).evaluate(q);

            // atom-water
            double aw_sum = std::inner_product(distance_profiles.aw.begin(ff1), distance_profiles.aw.end(ff1), sinqd_table->begin(q), 0.0);
            Iq[q-q0] += 2*free_params.cw*aw_sum*ff_table.index(ff1, form_factor::water_bin).evaluate(q);
        }

        // exv-exv
        double xx_sum = std::inner_product(distance_profiles.aa.begin(form_factor::exv_bin, form_factor::exv_bin), distance_profiles.aa.end(form_factor::exv_bin, form_factor::exv_bin), sinqd_table->begin(q), 0.0);
        Iq[q-q0] += cx*cx*xx_sum*ff_table.index(form_factor::exv_bin, form_factor::exv_bin).evaluate(q);

        // exv-water
        double ew_sum = std::inner_product(distance_profiles.aw.begin(form_factor::exv_bin), distance_profiles.aw.end(form_factor::exv_bin), sinqd_table->begin(q), 0.0);
        Iq[q-q0] -= 2*cx*free_params.cw*ew_sum*ff_table.index(form_factor::exv_bin, form_factor::water_bin).evaluate(q);

        // water-water
        double ww_sum = std::inner_product(distance_profiles.ww.begin(), distance_profiles.ww.end(), sinqd_table->begin(q), 0.0);
        Iq[q-q0] += free_params.cw*free_params.cw*ww_sum*ff_table.index(form_factor::water_bin, form_factor::water_bin).evaluate(q);
    }
    return ScatteringProfile(Iq, debye_axis);
}

template<typename FormFactorTableType>
SimpleDataset CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::debye_transform(const std::vector<double>&) const {
    throw except::not_implemented("CompositeDistanceHistogramFFGrid::debye_transform(const std::vector<double>& q) const");
}

template<typename FormFactorTableType>
const std::vector<double>& CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::get_counts() const {
    p = std::vector<double>(DistanceHistogram::get_counts().size(), 0);
    auto& p_pp = get_aa_counts();
    auto& p_hp = get_aw_counts();
    auto& p_hh = get_ww_counts();
    for (unsigned int i = 0; i < p.size(); ++i) {
        p[i] = p_pp.index(i) + 2*free_params.cw*p_hp.index(i) + free_params.cw*free_params.cw*p_hh.index(i);
    }
    return p.data;
}

template<typename FormFactorTableType>
Distribution1D& CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::get_aa_counts() {
    return const_cast<Distribution1D&>(const_cast<const CompositeDistanceHistogramFFAvgBase*>(this)->get_aa_counts());
}

template<typename FormFactorTableType>
const Distribution1D& CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::get_aa_counts() const {
    if (!cached_profiles.p_aa.empty()) {return cached_profiles.p_aa;}
    cached_profiles.p_aa = Distribution1D(axis.bins, 0);
    for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
        for (unsigned int ff2 = 0; ff2 < form_factor::get_count_without_excluded_volume(); ++ff2) {
            std::transform(cached_profiles.p_aa.begin(), cached_profiles.p_aa.end(), distance_profiles.aa.begin(ff1, ff2), cached_profiles.p_aa.begin(), std::plus<>());
        }
    }
    return cached_profiles.p_aa;
}

template<typename FormFactorTableType>
Distribution1D& CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::get_aw_counts() {return const_cast<Distribution1D&>(const_cast<const CompositeDistanceHistogramFFAvgBase*>(this)->get_aw_counts());}

template<typename FormFactorTableType>
const Distribution1D& CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::get_aw_counts() const {
    if (!cached_profiles.p_aw.empty()) {return cached_profiles.p_aw;}
    cached_profiles.p_aw = Distribution1D(axis.bins, 0);
    for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
        std::transform(cached_profiles.p_aw.begin(), cached_profiles.p_aw.end(), distance_profiles.aw.begin(ff1), cached_profiles.p_aw.begin(), std::plus<>());
    }
    return cached_profiles.p_aw;
}

template<typename FormFactorTableType>
Distribution1D& CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::get_ww_counts() {return const_cast<Distribution1D&>(const_cast<const CompositeDistanceHistogramFFAvgBase*>(this)->get_ww_counts());}

template<typename FormFactorTableType>
const Distribution1D& CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::get_ww_counts() const {
    if (!cached_profiles.p_ww.empty()) {return cached_profiles.p_ww;}
    cached_profiles.p_ww = Distribution1D(axis.bins, 0);
    std::transform(cached_profiles.p_ww.begin(), cached_profiles.p_ww.end(), distance_profiles.ww.begin(), cached_profiles.p_ww.begin(), std::plus<>());
    return cached_profiles.p_ww;
}

template<typename FormFactorTableType>
const Distribution3D& CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::get_aa_counts_ff() const {
    return distance_profiles.aa;
}

template<typename FormFactorTableType>
Distribution3D& CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::get_aa_counts_ff() {
    return distance_profiles.aa;
}

template<typename FormFactorTableType>
const Distribution2D& CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::get_aw_counts_ff() const {
    return distance_profiles.aw;
}

template<typename FormFactorTableType>
Distribution2D& CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::get_aw_counts_ff() {
    return distance_profiles.aw;
}

template<typename FormFactorTableType>
const Distribution1D& CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::get_ww_counts_ff() const {
    return distance_profiles.ww;
}

template<typename FormFactorTableType>
Distribution1D& CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::get_ww_counts_ff() {
    return distance_profiles.ww;
}

template<typename FormFactorTableType>
void CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::apply_water_scaling_factor(double k) {
    free_params.cw = k;
}

template<typename FormFactorTableType>
void CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::apply_excluded_volume_scaling_factor(double k) {
    free_params.cx = k;
}

template<typename FormFactorTableType>
ScatteringProfile CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::get_profile_aa() const {
    const auto& ff_table = get_ff_table();
    auto sinqd_table = get_sinc_table();
    Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
    unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

    std::vector<double> Iq(debye_axis.bins, 0);
    for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
        for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
            for (unsigned int ff2 = 0; ff2 < form_factor::get_count_without_excluded_volume(); ++ff2) {
                double aa_sum = std::inner_product(distance_profiles.aa.begin(ff1, ff2), distance_profiles.aa.end(ff1, ff2), sinqd_table->begin(q), 0.0);
                Iq[q-q0] += aa_sum*ff_table.index(ff1, ff2).evaluate(q);
            }
        }
    }
    return ScatteringProfile(std::move(Iq), debye_axis);
}

template<typename FormFactorTableType>
ScatteringProfile CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::get_profile_ax() const {
    const auto& ff_table = get_ff_table();
    auto sinqd_table = get_sinc_table();
    Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
    unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

    std::vector<double> Iq(debye_axis.bins, 0);
    for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
        double cx = exv_factor(q);
        for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
            double ax_sum = std::inner_product(distance_profiles.aa.begin(ff1, form_factor::exv_bin), distance_profiles.aa.end(ff1, form_factor::exv_bin), sinqd_table->begin(q), 0.0);
            Iq[q-q0] += 2*cx*ax_sum*ff_table.index(ff1, form_factor::exv_bin).evaluate(q);
        }
    }
    return ScatteringProfile(std::move(Iq), debye_axis);
}

template<typename FormFactorTableType>
ScatteringProfile CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::get_profile_xx() const {
    const auto& ff_table = get_ff_table();
    auto sinqd_table = get_sinc_table();
    Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
    unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

    std::vector<double> Iq(debye_axis.bins, 0);
    for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
        double cx = exv_factor(q);
        double xx_sum = std::inner_product(distance_profiles.aa.begin(form_factor::exv_bin, form_factor::exv_bin), distance_profiles.aa.end(form_factor::exv_bin, form_factor::exv_bin), sinqd_table->begin(q), 0.0);
        Iq[q-q0] += cx*cx*xx_sum*ff_table.index(form_factor::exv_bin, form_factor::exv_bin).evaluate(q);
    }
    return ScatteringProfile(std::move(Iq), debye_axis);
}

template<typename FormFactorTableType>
ScatteringProfile CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::get_profile_wx() const {
    const auto& ff_table = get_ff_table();
    auto sinqd_table = get_sinc_table();
    Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
    unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

    std::vector<double> Iq(debye_axis.bins, 0);
    for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
        double cx = exv_factor(q);
        double ew_sum = std::inner_product(distance_profiles.aw.begin(form_factor::exv_bin), distance_profiles.aw.end(form_factor::exv_bin), sinqd_table->begin(q), 0.0);
        Iq[q-q0] += 2*cx*free_params.cw*ew_sum*ff_table.index(form_factor::exv_bin, form_factor::water_bin).evaluate(q);
    }
    return ScatteringProfile(std::move(Iq), debye_axis);
}

template<typename FormFactorTableType>
ScatteringProfile CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::get_profile_aw() const {
    const auto& ff_table = get_ff_table();
    auto sinqd_table = get_sinc_table();
    Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
    unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

    std::vector<double> Iq(debye_axis.bins, 0);
    for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
        for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
            double aw_sum = std::inner_product(distance_profiles.aw.begin(ff1), distance_profiles.aw.end(ff1), sinqd_table->begin(q), 0.0);
            Iq[q-q0] += 2*free_params.cw*aw_sum*ff_table.index(ff1, form_factor::water_bin).evaluate(q);
        }
    }
    return ScatteringProfile(std::move(Iq), debye_axis);
}

template<typename FormFactorTableType>
ScatteringProfile CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::get_profile_ww() const {
    const auto& ff_table = get_ff_table();
    auto sinqd_table = get_sinc_table();
    Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
    unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

    std::vector<double> Iq(debye_axis.bins, 0);
    for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
        double ww_sum = std::inner_product(distance_profiles.ww.begin(), distance_profiles.ww.end(), sinqd_table->begin(q), 0.0);
        Iq[q-q0] += free_params.cw*free_params.cw*ww_sum*ff_table.index(form_factor::water_bin, form_factor::water_bin).evaluate(q);
    }
    return ScatteringProfile(std::move(Iq), debye_axis);
}

template class hist::CompositeDistanceHistogramFFAvgBase<form_factor::storage::atomic::table_t>;