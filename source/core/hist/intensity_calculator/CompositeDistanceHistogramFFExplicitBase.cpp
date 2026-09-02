// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/intensity_calculator/CompositeDistanceHistogramFFExplicitBase.h>
#include <form_factor/lookup/FormFactorManager.h>
#include <hist/Histogram.h>
#include <table/ArrayDebyeTable.h>
#include <form_factor/lookup/FormFactorProduct.h>
#include <settings/HistogramSettings.h>
#include <utility/MultiThreading.h>

using namespace ausaxs;
using namespace ausaxs::hist;

template<typename AA, typename AXFormFactorTableType, typename XX>
CompositeDistanceHistogramFFExplicitBase<AA, AXFormFactorTableType, XX>::CompositeDistanceHistogramFFExplicitBase() = default;

template<typename AA, typename AXFormFactorTableType, typename XX>
CompositeDistanceHistogramFFExplicitBase<AA, AXFormFactorTableType, XX>::CompositeDistanceHistogramFFExplicitBase(const CompositeDistanceHistogramFFExplicitBase&) = default;

template<typename AA, typename AX, typename XX>
CompositeDistanceHistogramFFExplicitBase<AA, AX, XX>::CompositeDistanceHistogramFFExplicitBase(CompositeDistanceHistogramFFExplicitBase&&) noexcept = default;

template<typename AA, typename AX, typename XX>
CompositeDistanceHistogramFFExplicitBase<AA, AX, XX>& CompositeDistanceHistogramFFExplicitBase<AA, AX, XX>::operator=(CompositeDistanceHistogramFFExplicitBase&&) noexcept = default;

template<typename AA, typename AX, typename XX>
CompositeDistanceHistogramFFExplicitBase<AA, AX, XX>& CompositeDistanceHistogramFFExplicitBase<AA, AX, XX>::operator=(const CompositeDistanceHistogramFFExplicitBase&) = default;

template<typename AA, typename AXFormFactorTableType, typename XX>
CompositeDistanceHistogramFFExplicitBase<AA, AXFormFactorTableType, XX>::~CompositeDistanceHistogramFFExplicitBase() = default;

template<typename AA, typename AXFormFactorTableType, typename XX>
CompositeDistanceHistogramFFExplicitBase<AA, AXFormFactorTableType, XX>::CompositeDistanceHistogramFFExplicitBase(
    hist::Distribution3D&& p_aa, 
    hist::Distribution2D&& p_aw, 
    hist::Distribution1D&& p_ww,
    hist::Distribution1D&& p_tot
) : CompositeDistanceHistogramFFAvgBase<AA>(std::move(p_aa), std::move(p_aw), std::move(p_ww), std::move(p_tot)) {}

template<typename AA, typename AXFormFactorTableType, typename XX>
CompositeDistanceHistogramFFExplicitBase<AA, AXFormFactorTableType, XX>::CompositeDistanceHistogramFFExplicitBase(
    hist::Distribution3D&& p_aa, 
    hist::Distribution2D&& p_aw, 
    hist::Distribution1D&& p_ww, 
    hist::WeightedDistribution1D&& p_tot
) : CompositeDistanceHistogramFFAvgBase<AA>(std::move(p_aa), std::move(p_aw), std::move(p_ww), std::move(p_tot))  {}

template<typename AA, typename AXFormFactorTableType, typename XX>
const AA CompositeDistanceHistogramFFExplicitBase<AA, AXFormFactorTableType, XX>::get_ffaa_table() const {
    return this->get_ff_table();
}

template<typename AA, typename AXFormFactorTableType, typename XX>
double CompositeDistanceHistogramFFExplicitBase<AA, AXFormFactorTableType, XX>::exv_factor(double q, double cx) {
    constexpr double rm2 = constants::radius::average_atomic_radius*constants::radius::average_atomic_radius/4;
    return std::pow(cx, 3)*std::exp(-rm2*(std::pow(cx, 2) - 1)*q*q);
}

template<typename AA, typename AXFormFactorTableType, typename XX>
double CompositeDistanceHistogramFFExplicitBase<AA, AXFormFactorTableType, XX>::exv_factor(double q) const {
    return exv_factor(q, this->free_params.cx);
}

template<typename AA, typename AXFormFactorTableType, typename XX>
void CompositeDistanceHistogramFFExplicitBase<AA, AXFormFactorTableType, XX>::cache_refresh_intensity_exv(const std::vector<double>& cx, bool cw_changed, bool cx_changed) const {
    auto pool = utility::multi_threading::get_global_pool();

    unsigned int bins = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax).bins;
    unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin);

    // these lazily initialize shared state, so they must be resolved before any job is submitted
    const auto* ff_ax_table = &get_ffax_table();
    const auto* ff_xx_table = &get_ffxx_table();

    // aa, ax and xx all share the same distance histogram, so the atomic sinqd cache is reused with the
    // ax and xx form factor tables. Likewise wx reuses the aw cache.
    if (cx_changed) {
        // ax
        pool->detach_task([this, &cx, q0, bins, ff_ax_table] () {
            for (unsigned int ff1 = form_factor::start_index_for_explicit_exv(); ff1 < form_factor::get_active_count(); ++ff1) {
                for (unsigned int ff2 = form_factor::start_index_for_explicit_exv(); ff2 < form_factor::get_active_count(); ++ff2) {
                    for (unsigned int q = q0; q < q0+bins; ++q) {
                        this->cache.intensity_profiles.ax[q-q0] += 
                            this->free_params.crho*cx[q-q0]*this->cache.sinqd.aa.index(ff1, ff2, q-q0)
                            *(ff_ax_table->index(ff1, ff2).evaluate(q) + ff_ax_table->index(ff2, ff1).evaluate(q))
                        ;
                    }
                }
            }
        });

        // xx
        pool->detach_task([this, &cx, q0, bins, ff_xx_table] () {
            for (unsigned int ff1 = form_factor::start_index_for_explicit_exv(); ff1 < form_factor::get_active_count(); ++ff1) {
                for (unsigned int ff2 = form_factor::start_index_for_explicit_exv(); ff2 < form_factor::get_active_count(); ++ff2) {
                    for (unsigned int q = q0; q < q0+bins; ++q) {
                        this->cache.intensity_profiles.xx[q-q0] += 
                            std::pow(cx[q-q0]*this->free_params.crho, 2)*this->cache.sinqd.aa.index(ff1, ff2, q-q0)*ff_xx_table->index(ff1, ff2).evaluate(q)
                        ;
                    }
                }
            }
        });
    }

    if (cw_changed || cx_changed) {
        // wx
        pool->detach_task([this, &cx, q0, bins, ff_ax_table] () {
            for (unsigned int ff1 = form_factor::start_index_for_explicit_exv(); ff1 < form_factor::get_active_count(); ++ff1) {
                for (unsigned int q = q0; q < q0+bins; ++q) {
                    // only the atom carries excluded volume, so the exv form factor is evaluated only for the atom, not the water. 
                    // See the CompositeDistanceHistogramFFAvgBase class documentation.
                    this->cache.intensity_profiles.wx[q-q0] += 
                        2*this->free_params.crho*cx[q-q0]*this->free_params.cw*this->cache.sinqd.aw.index(ff1, q-q0)
                        *ff_ax_table->index(form_factor::water_bin, ff1).evaluate(q)
                    ;
                }
            }
        });
    }
}

template class hist::CompositeDistanceHistogramFFExplicitBase<
    form_factor::lookup::table_t, form_factor::lookup::table_t, form_factor::lookup::table_t
>;