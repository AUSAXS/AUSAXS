// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/intensity_calculator/CompositeDistanceHistogramFFGrid.h>
#include <form_factor/FormFactorType.h>
#include <form_factor/lookup/FormFactorManager.h>
#include <utility/MultiThreading.h>
#include <settings/HistogramSettings.h>

using namespace ausaxs;
using namespace ausaxs::hist;
using namespace ausaxs::form_factor;

CompositeDistanceHistogramFFGrid::CompositeDistanceHistogramFFGrid(CompositeDistanceHistogramFFGrid&&) noexcept = default;
CompositeDistanceHistogramFFGrid& CompositeDistanceHistogramFFGrid::operator=(CompositeDistanceHistogramFFGrid&&) noexcept = default;
CompositeDistanceHistogramFFGrid::~CompositeDistanceHistogramFFGrid() = default;

CompositeDistanceHistogramFFGrid::CompositeDistanceHistogramFFGrid(
    hist::Distribution3D&& p_aa, 
    hist::Distribution2D&& p_aw, 
    hist::Distribution1D&& p_ww, 
    hist::WeightedDistribution1D&& p_tot_aa,
    hist::WeightedDistribution1D&& p_tot_ax,
    hist::WeightedDistribution1D&& p_tot_xx
) : CompositeDistanceHistogramFFGridBase(std::move(p_aa), std::move(p_aw), std::move(p_ww), std::move(p_tot_aa)) {
    ff_table = generate_ff_table();
    initialize_grid_axes(p_tot_ax.get_weighted_axis(), p_tot_xx.get_weighted_axis());
}

double CompositeDistanceHistogramFFGrid::exv_factor(double) const {
    return free_params.cx;
}

void CompositeDistanceHistogramFFGrid::cache_refresh_sinqd_exv() const {
    auto pool = utility::multi_threading::get_global_pool();
    const auto& sinqd_table_ax = get_sinc_table_ax();
    const auto& sinqd_table_xx = get_sinc_table_xx();

    Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
    unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin);

    if (cache.sinqd.ax.empty()) {
        cache.sinqd.ax = container::Container2D<double>(form_factor::get_total_ff_count(), debye_axis.bins);
        cache.sinqd.xx = container::Container1D<double>(debye_axis.bins);
        cache.sinqd.wx = container::Container1D<double>(debye_axis.bins);
    }

    for (unsigned int ff1 = form_factor::start_index_for_explicit_exv(); ff1 < form_factor::get_active_count(); ++ff1) {
        pool->detach_task([this, q0, bins=debye_axis.bins, ff1, sinqd_table_ax] () {
            for (unsigned int q = q0; q < q0+bins; ++q) {
                cache.sinqd.ax.index(ff1, q-q0) = 2*std::inner_product(distance_profiles.aa.begin(ff1, form_factor::exv_bin), distance_profiles.aa.end(ff1, form_factor::exv_bin), sinqd_table_ax->begin(q), 0.0);
            }
        });
    }
    pool->detach_task([&] () {
        for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
            cache.sinqd.xx.index(q-q0) = std::inner_product(distance_profiles.aa.begin(form_factor::exv_bin, form_factor::exv_bin), distance_profiles.aa.end(form_factor::exv_bin, form_factor::exv_bin), sinqd_table_xx->begin(q), 0.0);
            cache.sinqd.wx.index(q-q0) = 2*std::inner_product(distance_profiles.aw.begin(form_factor::exv_bin), distance_profiles.aw.end(form_factor::exv_bin), sinqd_table_ax->begin(q), 0.0);
        }
    });
    pool->wait();
}
