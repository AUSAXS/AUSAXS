// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/histogram_manager/detail/GridExvHelpers.h>

#include <form_factor/FormFactorType.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFAvg.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFGrid.h>

#include <algorithm>

using namespace ausaxs;
using namespace ausaxs::hist;
using namespace ausaxs::hist::detail::grid_exv;

AtomicDistributions AtomicDistributions::take(CompositeDistanceHistogramFFAvg& base) {
    // the total is copied before anything is moved out, since it is derived from the rest
    WeightedDistribution1D p_tot = base.get_weighted_counts();
    p_tot.set_bin_centers(base.get_d_axis());
    return AtomicDistributions{
        .p_aa = std::move(base.get_raw_aa_counts_by_ff()),
        .p_aw = std::move(base.get_raw_aw_counts_by_ff()),
        .p_ww = std::move(base.get_raw_ww_counts_by_ff()),
        .p_tot = std::move(p_tot)
    };
}

int AtomicDistributions::grow(int bins) {
    if (bins <= p_tot.size()) {return p_tot.size();}
    p_aa.resize(bins);
    p_aw.resize(bins);
    p_ww.resize(bins);
    p_tot.resize(bins);
    return bins;
}

std::unique_ptr<CompositeDistanceHistogramFFGrid> hist::detail::grid_exv::splice(
    AtomicDistributions atomic, const WeightedDistribution2D& p_ax, const WeightedDistribution1D& p_wx, WeightedDistribution1D p_xx
) {
    int bins = atomic.p_tot.size();

    // calculate weighted distance bins
    WeightedDistribution1D p_tot_ax(std::max<int>(bins, p_wx.size()));
    for (int i = 0; i < bins; ++i) {
        p_tot_ax.add_index(i, p_wx.index(i));
    }
    for (int i = 0; i < p_ax.size_x(); ++i) {
        for (int j = 0; j < bins; ++j) {
            p_tot_ax.add_index(j, p_ax.index(i, j));
        }
    }

    // overwrite the excluded volume of the averaged model with the grid-based one
    Distribution2D ax(p_ax);
    Distribution1D wx(p_wx);
    Distribution1D xx(p_xx);
    for (int i = 0; i < atomic.p_aa.size_x(); ++i) {
        std::move(ax.begin(i), ax.begin(i)+bins, atomic.p_aa.begin(i, form_factor::exv_bin));
    }
    std::move(wx.begin(), wx.begin()+bins, atomic.p_aw.begin(form_factor::exv_bin));
    std::move(xx.begin(), xx.begin()+bins, atomic.p_aa.begin(form_factor::exv_bin, form_factor::exv_bin));

    return std::make_unique<CompositeDistanceHistogramFFGrid>(
        std::move(atomic.p_aa),
        std::move(atomic.p_aw),
        std::move(atomic.p_ww),
        std::move(atomic.p_tot),
        std::move(p_tot_ax),
        std::move(p_xx)
    );
}
