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

void CompositeDistanceHistogramFFAvg::cache_refresh_sinqd_exv() const {
    Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
    unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin);

    if (cache.sinqd.ax.empty()) {
        cache.sinqd.ax = container::Container2D<double>(form_factor::get_total_ff_count(), debye_axis.bins);
        cache.sinqd.xx = container::Container1D<double>(debye_axis.bins);
        cache.sinqd.wx = container::Container1D<double>(debye_axis.bins);
    }

    const unsigned int ff_start = form_factor::start_index_for_explicit_exv();
    const unsigned int n_active = form_factor::get_active_count();
    const double Z = Z_exv_avg;

    for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
        double total = 0;
        for (unsigned int a = ff_start; a < n_active; ++a) {
            double row = 0, col = 0;
            for (unsigned int b = ff_start; b < n_active; ++b) {
                row += cache.sinqd.aa.index(a, b, q-q0);
                col += cache.sinqd.aa.index(b, a, q-q0);
            }
            cache.sinqd.ax.index(a, q-q0) = Z*(row + col);
            total += row;
        }
        cache.sinqd.xx.index(q-q0) = Z*Z*total;

        double wx = 0;
        for (unsigned int a = ff_start; a < n_active; ++a) {
            wx += cache.sinqd.aw.index(a, q-q0);
        }
        cache.sinqd.wx.index(q-q0) = 2*Z*wx;
    }
}
