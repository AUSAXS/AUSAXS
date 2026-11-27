// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/histogram_manager/detail/GridExvHelpers.h>

#include <form_factor/FormFactorType.h>
#include <hist/scattering_model/ScatteringModelGrid.h>

#include <algorithm>

using namespace ausaxs;
using namespace ausaxs::hist;
using namespace ausaxs::hist::detail::grid_exv;

AtomicDistributions AtomicDistributions::from(ManagerDistributions<true, true>&& d) { // NOLINT
    return AtomicDistributions{
        .p_aa = Distribution3D<hist::Shape::Triangular>(d.p_aa),
        .p_aw = Distribution2D(d.p_aw),
        .p_ww = Distribution1D(d.p_ww),
        .self_aa = std::move(d.self_aa),
        .self_ww = d.self_ww,
        .p_tot = std::move(d.p_tot)
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

std::unique_ptr<ScatteringModelGrid> hist::detail::grid_exv::splice(
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

    // write the grid-based excluded volume into the excluded volume slot
    Distribution2D ax(p_ax);
    Distribution1D wx(p_wx);
    Distribution1D xx(p_xx);
    for (int i = 0; i < atomic.p_aa.size_x(); ++i) {
        std::move(ax.begin(i), ax.begin(i)+bins, atomic.p_aa.begin(i, form_factor::exv_bin));
    }
    std::move(wx.begin(), wx.begin()+bins, atomic.p_aw.begin(form_factor::exv_bin));
    std::move(xx.begin(), xx.begin()+bins, atomic.p_aa.begin(form_factor::exv_bin, form_factor::exv_bin));

    return std::make_unique<ScatteringModelGrid>(
        std::move(atomic.p_aa),
        std::move(atomic.p_aw),
        std::move(atomic.p_ww),
        std::move(atomic.self_aa),
        atomic.self_ww,
        std::move(atomic.p_tot),
        std::move(p_tot_ax),
        std::move(p_xx)
    );
}
