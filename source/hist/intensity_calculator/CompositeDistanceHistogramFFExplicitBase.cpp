/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <hist/intensity_calculator/CompositeDistanceHistogramFFExplicitBase.h>
#include <hist/Histogram.h>
#include <table/ArrayDebyeTable.h>
#include <form_factor/FormFactor.h>
#include <form_factor/PrecalculatedFormFactorProduct.h>
#include <form_factor/PrecalculatedExvFormFactorProduct.h>
#include <settings/HistogramSettings.h>
#include <constants/ConstantsMath.h>
#include <math/ConstexprMath.h>

using namespace hist;

template<typename AAFormFactorTableType, typename AXFormFactorTableType, typename XXFormFactorTableType>
CompositeDistanceHistogramFFExplicitBase<AAFormFactorTableType, AXFormFactorTableType, XXFormFactorTableType>::CompositeDistanceHistogramFFExplicitBase() = default;

template<typename AAFormFactorTableType, typename AXFormFactorTableType, typename XXFormFactorTableType>
CompositeDistanceHistogramFFExplicitBase<AAFormFactorTableType, AXFormFactorTableType, XXFormFactorTableType>::CompositeDistanceHistogramFFExplicitBase(
    hist::Distribution3D&& p_aa, 
    hist::Distribution3D&& p_ax, 
    hist::Distribution3D&& p_xx, 
    hist::Distribution2D&& p_aw, 
    hist::Distribution2D&& p_wx, 
    hist::Distribution1D&& p_ww,
    hist::Distribution1D&& p_tot
) : CompositeDistanceHistogramFFAvgBase<AAFormFactorTableType>(std::move(p_aa), std::move(p_aw), std::move(p_ww), std::move(p_tot)), cp_ax(std::move(p_ax)), cp_xx(std::move(p_xx)), cp_wx(std::move(p_wx)) {}

template<typename AAFormFactorTableType, typename AXFormFactorTableType, typename XXFormFactorTableType>
CompositeDistanceHistogramFFExplicitBase<AAFormFactorTableType, AXFormFactorTableType, XXFormFactorTableType>::CompositeDistanceHistogramFFExplicitBase(
    hist::Distribution3D&& p_aa, 
    hist::Distribution3D&& p_ax, 
    hist::Distribution3D&& p_xx, 
    hist::Distribution2D&& p_aw, 
    hist::Distribution2D&& p_wx, 
    hist::Distribution1D&& p_ww, 
    hist::WeightedDistribution1D&& p_tot
) : CompositeDistanceHistogramFFAvgBase<AAFormFactorTableType>(std::move(p_aa), std::move(p_aw), std::move(p_ww), std::move(p_tot)), cp_ax(std::move(p_ax)), cp_xx(std::move(p_xx)), cp_wx(std::move(p_wx)) {}

template<typename AAFormFactorTableType, typename AXFormFactorTableType, typename XXFormFactorTableType>
CompositeDistanceHistogramFFExplicitBase<AAFormFactorTableType, AXFormFactorTableType, XXFormFactorTableType>::~CompositeDistanceHistogramFFExplicitBase() = default;


template<typename AAFormFactorTableType, typename AXFormFactorTableType, typename XXFormFactorTableType>
const AAFormFactorTableType CompositeDistanceHistogramFFExplicitBase<AAFormFactorTableType, AXFormFactorTableType, XXFormFactorTableType>::get_ffaa_table() const {
    return this->get_ff_table();
}

template<typename AAFormFactorTableType, typename AXFormFactorTableType, typename XXFormFactorTableType>
double CompositeDistanceHistogramFFExplicitBase<AAFormFactorTableType, AXFormFactorTableType, XXFormFactorTableType>::exv_factor(double q) const {
    constexpr double rm2 = constants::form_factor::sigma_excluded_volume*constants::form_factor::sigma_excluded_volume;
    return std::pow(this->cx, 3)*std::exp(-rm2*(this->cx*this->cx - 1)*q*q/4);
}

template<typename AAFormFactorTableType, typename AXFormFactorTableType, typename XXFormFactorTableType>
ScatteringProfile CompositeDistanceHistogramFFExplicitBase<AAFormFactorTableType, AXFormFactorTableType, XXFormFactorTableType>::debye_transform() const {
    const auto& ff_aa_table = get_ffaa_table();
    const auto& ff_ax_table = get_ffax_table();
    const auto& ff_xx_table = get_ffxx_table();
    auto sinqd_table = this->get_sinc_table();

    // calculate the Debye scattering intensity
    Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
    unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

    std::vector<double> Iq(debye_axis.bins, 0);
    for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
        double cx = exv_factor(constants::axes::q_vals[q]);
        for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
            for (unsigned int ff2 = 0; ff2 < form_factor::get_count_without_excluded_volume(); ++ff2) {
                // atom-atom
                double aa_sum = std::inner_product(this->cp_aa.begin(ff1, ff2), this->cp_aa.end(ff1, ff2), sinqd_table->begin(q), 0.0);
                Iq[q-q0] += aa_sum*ff_aa_table.index(ff1, ff2).evaluate(q);

                // atom-exv
                double ax_sum = std::inner_product(cp_ax.begin(ff1, ff2), cp_ax.end(ff1, ff2), sinqd_table->begin(q), 0.0);
                Iq[q-q0] -= 2*cx*ax_sum*ff_ax_table.index(ff1, ff2).evaluate(q);

                // exv-exv
                double xx_sum = std::inner_product(cp_xx.begin(ff1, ff2), cp_xx.end(ff1, ff2), sinqd_table->begin(q), 0.0);
                Iq[q-q0] += cx*cx*xx_sum*ff_xx_table.index(ff1, ff2).evaluate(q);
            }

            // atom-water
            double aw_sum = std::inner_product(this->cp_aw.begin(ff1), this->cp_aw.end(ff1), sinqd_table->begin(q), 0.0);
            Iq[q-q0] += 2*this->cw*aw_sum*ff_aa_table.index(ff1, form_factor::water_bin).evaluate(q);

            // exv-water
            double wx_sum = std::inner_product(cp_wx.begin(ff1), cp_wx.end(ff1), sinqd_table->begin(q), 0.0);
            Iq[q-q0] -= 2*cx*this->cw*wx_sum*ff_ax_table.index(form_factor::water_bin, ff1).evaluate(q);
        }

        // water-water
        double ww_sum = std::inner_product(this->cp_ww.begin(), this->cp_ww.end(), sinqd_table->begin(q), 0.0);
        Iq[q-q0] += this->cw*this->cw*ww_sum*ff_aa_table.index(form_factor::water_bin, form_factor::water_bin).evaluate(q);
    }
    return ScatteringProfile(std::move(Iq), debye_axis);
}

template<typename AAFormFactorTableType, typename AXFormFactorTableType, typename XXFormFactorTableType>
SimpleDataset CompositeDistanceHistogramFFExplicitBase<AAFormFactorTableType, AXFormFactorTableType, XXFormFactorTableType>::debye_transform(const std::vector<double>&) const {
    throw except::not_implemented("CompositeDistanceHistogramFFGrid::debye_transform(const std::vector<double>& q) const");
}

template<typename AAFormFactorTableType, typename AXFormFactorTableType, typename XXFormFactorTableType>
ScatteringProfile CompositeDistanceHistogramFFExplicitBase<AAFormFactorTableType, AXFormFactorTableType, XXFormFactorTableType>::get_profile_ax() const {
    const auto& ff_ax_table = get_ffax_table();
    auto sinqd_table = this->get_sinc_table();
    Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
    unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

    std::vector<double> Iq(debye_axis.bins, 0);
    for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
        double cx = exv_factor(constants::axes::q_vals[q]);
        for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
            for (unsigned int ff2 = 0; ff2 < form_factor::get_count_without_excluded_volume(); ++ff2) {
                double ax_sum = std::inner_product(cp_ax.begin(ff1, ff2), cp_ax.end(ff1, ff2), sinqd_table->begin(q), 0.0);
                Iq[q-q0] += 2*cx*ax_sum*ff_ax_table.index(ff1, ff2).evaluate(q);
            }
        }
    }
    return ScatteringProfile(std::move(Iq), debye_axis);
}

template<typename AAFormFactorTableType, typename AXFormFactorTableType, typename XXFormFactorTableType>
ScatteringProfile CompositeDistanceHistogramFFExplicitBase<AAFormFactorTableType, AXFormFactorTableType, XXFormFactorTableType>::get_profile_xx() const {
    const auto& ff_xx_table = get_ffxx_table();
    auto sinqd_table = this->get_sinc_table();
    Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
    unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

    std::vector<double> Iq(debye_axis.bins, 0);
    for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
        double cx2 = std::pow(exv_factor(constants::axes::q_vals[q]), 2);
        for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
            for (unsigned int ff2 = 0; ff2 < form_factor::get_count_without_excluded_volume(); ++ff2) {
                double xx_sum = std::inner_product(cp_xx.begin(ff1, ff2), cp_xx.end(ff1, ff2), sinqd_table->begin(q), 0.0);
                Iq[q-q0] += cx2*xx_sum*ff_xx_table.index(ff1, ff2).evaluate(q);
            }
        }
    }
    return ScatteringProfile(std::move(Iq), debye_axis);
}

template<typename AAFormFactorTableType, typename AXFormFactorTableType, typename XXFormFactorTableType>
ScatteringProfile CompositeDistanceHistogramFFExplicitBase<AAFormFactorTableType, AXFormFactorTableType, XXFormFactorTableType>::get_profile_wx() const {
    const auto& ff_ax_table = get_ffax_table();
    auto sinqd_table = this->get_sinc_table();
    Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
    unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

    std::vector<double> Iq(debye_axis.bins, 0);
    unsigned int ff_w_index = static_cast<int>(form_factor::form_factor_t::OH);
    for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
        double cx = exv_factor(constants::axes::q_vals[q]);
        for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
            double wx_sum = std::inner_product(cp_wx.begin(ff1), cp_wx.end(ff1), sinqd_table->begin(q), 0.0);
            Iq[q-q0] += 2*cx*this->cw*wx_sum*ff_ax_table.index(ff_w_index, ff1).evaluate(q);
        }
    }
    return ScatteringProfile(std::move(Iq), debye_axis);
}

template class hist::CompositeDistanceHistogramFFExplicitBase<form_factor::storage::atomic::table_t, form_factor::storage::cross::table_t, form_factor::storage::exv::table_t>;