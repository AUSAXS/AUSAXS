// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/intensity_calculator/CompositeDistanceHistogramFFAvg.h>
#include <hist/distribution/Distribution1D.h>
#include <hist/distribution/Distribution2D.h>
#include <hist/distribution/Distribution3D.h>
#include <hist/distribution/WeightedDistribution1D.h>
#include <form_factor/FormFactorType.h>
#include <form_factor/lookup/FormFactorManager.h>
#include <settings/HistogramSettings.h>
#include <utility/MultiThreading.h>

using namespace ausaxs;
using namespace ausaxs::hist;

CompositeDistanceHistogramFFAvg::CompositeDistanceHistogramFFAvg(
    hist::Distribution3D&& p_aa,
    hist::Distribution2D&& p_aw,
    hist::Distribution1D&& p_ww,
    hist::Distribution1D&& p_tot,
    double Z_exv_avg
) : CompositeDistanceHistogramFFAvgBase(std::move(p_aa), std::move(p_aw), std::move(p_ww), std::move(p_tot)), Z_exv_avg(Z_exv_avg) {}

CompositeDistanceHistogramFFAvg::CompositeDistanceHistogramFFAvg(
    hist::Distribution3D&& p_aa,
    hist::Distribution2D&& p_aw,
    hist::Distribution1D&& p_ww,
    hist::WeightedDistribution1D&& p_tot,
    double Z_exv_avg
) : CompositeDistanceHistogramFFAvgBase(std::move(p_aa), std::move(p_aw), std::move(p_ww), std::move(p_tot)), Z_exv_avg(Z_exv_avg) {}

void CompositeDistanceHistogramFFAvg::cache_refresh_intensity_exv(const std::vector<double>& cx, bool cw_changed, bool cx_changed) const {
    auto pool = utility::multi_threading::get_global_pool();

    unsigned int bins = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax).bins;
    unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin);
    const double Z = Z_exv_avg;

    // this lazily initializes shared state, so it must be resolved before any job is submitted
    const auto* ff_table = &get_ff_table();

    if (cx_changed) {
        // ax
        pool->detach_task([this, &cx, q0, bins, Z, ff_table] () {
            for (unsigned int ff1 = form_factor::start_index_for_explicit_exv(); ff1 < form_factor::get_active_count(); ++ff1) {
                for (unsigned int ff2 = form_factor::start_index_for_explicit_exv(); ff2 < form_factor::get_active_count(); ++ff2) {
                    for (unsigned int q = q0; q < q0+bins; ++q) {
                        cache.intensity_profiles.ax[q-q0] += Z*free_params.crho*cx[q-q0]*cache.sinqd.aa.index(ff1, ff2, q-q0)*(
                            ff_table->index(ff1, form_factor::exv_bin).evaluate(q) + ff_table->index(ff2, form_factor::exv_bin).evaluate(q)
                        );
                    }
                }
            }
        });

        // xx
        pool->detach_task([this, &cx, q0, bins, Z, ff_table] () {
            for (unsigned int ff1 = form_factor::start_index_for_explicit_exv(); ff1 < form_factor::get_active_count(); ++ff1) {
                for (unsigned int ff2 = form_factor::start_index_for_explicit_exv(); ff2 < form_factor::get_active_count(); ++ff2) {
                    for (unsigned int q = q0; q < q0+bins; ++q) {
                        cache.intensity_profiles.xx[q-q0] += std::pow(Z*cx[q-q0]*free_params.crho, 2)*cache.sinqd.aa.index(ff1, ff2, q-q0)
                            *ff_table->index(form_factor::exv_bin, form_factor::exv_bin).evaluate(q);
                    }
                }
            }
        });
    }

    if (cw_changed || cx_changed) {
        // wx
        pool->detach_task([this, &cx, q0, bins, Z, ff_table] () {
            for (unsigned int ff1 = form_factor::start_index_for_explicit_exv(); ff1 < form_factor::get_active_count(); ++ff1) {
                for (unsigned int q = q0; q < q0+bins; ++q) {
                    cache.intensity_profiles.wx[q-q0] += 2*Z*free_params.crho*cx[q-q0]*free_params.cw*cache.sinqd.aw.index(ff1, q-q0)
                        *ff_table->index(form_factor::exv_bin, form_factor::water_bin).evaluate(q);
                }
            }
        });
    }
}
