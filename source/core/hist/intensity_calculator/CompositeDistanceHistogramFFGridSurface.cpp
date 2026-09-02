// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/intensity_calculator/CompositeDistanceHistogramFFGridSurface.h>
#include <form_factor/lookup/FormFactorProduct.h>
#include <form_factor/lookup/FormFactorManager.h>
#include <form_factor/ExvFormFactor.h>
#include <table/ArrayDebyeTable.h>
#include <settings/GridSettings.h>
#include <settings/HistogramSettings.h>
#include <utility/MultiThreading.h>
#include <utility/Exceptions.h>
#include <dataset/SimpleDataset.h>

using namespace ausaxs;
using namespace ausaxs::hist;
using namespace ausaxs::form_factor;

CompositeDistanceHistogramFFGridSurface::CompositeDistanceHistogramFFGridSurface(CompositeDistanceHistogramFFGridSurface&&) noexcept = default;
CompositeDistanceHistogramFFGridSurface& CompositeDistanceHistogramFFGridSurface::operator=(CompositeDistanceHistogramFFGridSurface&&) noexcept = default;
CompositeDistanceHistogramFFGridSurface::~CompositeDistanceHistogramFFGridSurface() = default;

CompositeDistanceHistogramFFGridSurface::XXContainer CompositeDistanceHistogramFFGridSurface::XXContainer::operator+=(const CompositeDistanceHistogramFFGridSurface::XXContainer& other) {
    std::transform(interior.begin(), interior.end(), other.interior.begin(), interior.begin(), std::plus<>());
    std::transform(surface.begin(), surface.end(), other.surface.begin(), surface.begin(), std::plus<>());
    std::transform(cross.begin(), cross.end(), other.cross.begin(), cross.begin(), std::plus<>());
    return *this;
}

CompositeDistanceHistogramFFGridSurface::AXContainer CompositeDistanceHistogramFFGridSurface::AXContainer::operator+=(const CompositeDistanceHistogramFFGridSurface::AXContainer& other) {
    std::transform(interior.begin(), interior.end(), other.interior.begin(), interior.begin(), [](auto& a, const auto& b) {return a+b;});
    std::transform(surface.begin(), surface.end(), other.surface.begin(), surface.begin(), [](auto& a, const auto& b) {return a+b;});
    return *this;
}

CompositeDistanceHistogramFFGridSurface::WXContainer CompositeDistanceHistogramFFGridSurface::WXContainer::operator+=(const CompositeDistanceHistogramFFGridSurface::WXContainer& other) {
    std::transform(interior.begin(), interior.end(), other.interior.begin(), interior.begin(), std::plus<>());
    std::transform(surface.begin(), surface.end(), other.surface.begin(), surface.begin(), std::plus<>());
    return *this;
}

CompositeDistanceHistogramFFGridSurface::CompositeDistanceHistogramFFGridSurface(
    hist::Distribution3D&& p_aa, 
    hist::Distribution2D&& p_aw, 
    hist::Distribution1D&& p_ww, 
    XXContainer&& xx,
    AXContainer&& ax,
    WXContainer&& wx,
    hist::WeightedDistribution1D&& p_tot_aa,
    hist::WeightedDistribution1D&& p_tot_ax,
    hist::WeightedDistribution1D&& p_tot_xx
) : hist::CompositeDistanceHistogramFFGridBase(std::move(p_aa), std::move(p_aw), std::move(p_ww), std::move(p_tot_aa)),
    exv_distance_profiles{hist::Distribution1D(std::move(xx.interior)), hist::Distribution1D(std::move(xx.surface)), hist::Distribution1D(std::move(xx.cross)), 
                          hist::Distribution1D(std::move(wx.interior)), hist::Distribution1D(std::move(wx.surface)), hist::Distribution2D(std::move(ax.interior)), 
                          hist::Distribution2D(std::move(ax.surface))}
{
    initialize(p_tot_ax.get_weighted_axis(), p_tot_xx.get_weighted_axis());
}

Limit CompositeDistanceHistogramFFGridSurface::get_excluded_volume_scaling_factor_limits() const {
    return {0, 2};
}

double CompositeDistanceHistogramFFGridSurface::exv_factor(double q, double cx) {
    constexpr double rm2 = constants::radius::average_atomic_radius*constants::radius::average_atomic_radius/4;
    return std::pow(cx, 3)*std::exp(-rm2*(std::pow(cx, 2) - 1)*q*q);
}

double CompositeDistanceHistogramFFGridSurface::exv_factor(double q) const {
    return exv_factor(q, free_params.cx);
}

hist::Distribution1D CompositeDistanceHistogramFFGridSurface::evaluate_xx_distance_profile(double cx) const {
    hist::Distribution1D xx = exv_distance_profiles.xx_i;
    std::transform(xx.begin(), xx.end(), exv_distance_profiles.xx_s.begin(), xx.begin(), [cx] (double a, double b) {return a + std::pow(cx, 2)*b;});
    std::transform(xx.begin(), xx.end(), exv_distance_profiles.xx_c.begin(), xx.begin(), [cx] (double a, double b) {return a + cx*b;});
    return xx;
}

hist::Distribution1D CompositeDistanceHistogramFFGridSurface::evaluate_wx_distance_profile(double cx) const {
    hist::Distribution1D wx = exv_distance_profiles.wx_i;
    std::transform(wx.begin(), wx.end(), exv_distance_profiles.wx_s.begin(), wx.begin(), [cx] (double a, double b) {return a + cx*b;});
    return wx;
}

hist::Distribution2D CompositeDistanceHistogramFFGridSurface::evaluate_ax_distance_profile(double cx) const {
    hist::Distribution2D ax = exv_distance_profiles.ax_i;
    std::transform(ax.begin(), ax.end(), exv_distance_profiles.ax_s.begin(), ax.begin(), [cx] (double a, double b) {return a + cx*b;});
    return ax;
}

void CompositeDistanceHistogramFFGridSurface::initialize(std::vector<double>&& d_axis_ax, std::vector<double>&& d_axis_xx) {
    regenerate_ff_table();
    initialize_grid_axes(std::move(d_axis_ax), std::move(d_axis_xx));

    // fix the aa counts to also contain the exv contributions
    auto xx = evaluate_xx_distance_profile(1);
    auto wx = evaluate_wx_distance_profile(1);
    auto ax = evaluate_ax_distance_profile(1);

    auto& aa = CompositeDistanceHistogramFFAvgBase::get_raw_aa_counts_by_ff();
    for (unsigned int ff = form_factor::start_index_for_explicit_exv(); ff < form_factor::get_active_count(); ++ff) {
        std::transform(aa.begin(ff, form_factor::exv_bin), aa.end(ff, form_factor::exv_bin), ax.begin(ff), aa.begin(ff, form_factor::exv_bin), std::plus<double>());
    }
    std::transform(aa.begin(form_factor::exv_bin, form_factor::exv_bin), aa.end(form_factor::exv_bin, form_factor::exv_bin), xx.begin(), aa.begin(form_factor::exv_bin, form_factor::exv_bin), std::plus<double>());

    auto& aw = CompositeDistanceHistogramFFAvgBase::get_raw_aw_counts_by_ff();
    std::transform(aw.begin(form_factor::exv_bin), aw.end(form_factor::exv_bin), wx.begin(), aw.begin(form_factor::exv_bin), std::plus<double>());
}

void CompositeDistanceHistogramFFGridSurface::cache_refresh_intensity_exv(const std::vector<double>& cx, bool cw_changed, bool cx_changed) const {
    auto pool = utility::multi_threading::get_global_pool();

    unsigned int bins = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax).bins;
    unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin);

    // these lazily initialize shared state, so they must be resolved before any job is submitted
    const auto* ff_table = &get_ff_table();
    auto sinqd_table_ax = get_sinc_table_ax();
    auto sinqd_table_xx = get_sinc_table_xx();

    if (cx_changed) {
        // ax
        pool->detach_task([this, &cx, q0, bins, ff_table, sinqd_table_ax] () {
            for (unsigned int q = q0; q < q0+bins; ++q) {
                auto ax = evaluate_ax_distance_profile(cx[q-q0]);
                for (unsigned int ff1 = form_factor::start_index_for_explicit_exv(); ff1 < form_factor::get_active_count(); ++ff1) {
                    double ax_sum = std::inner_product(ax.begin(ff1), ax.end(ff1), sinqd_table_ax->begin(q), 0.0);
                    cache.intensity_profiles.ax[q-q0] += 2*free_params.crho*ax_sum*ff_table->index(ff1, form_factor::exv_bin).evaluate(q);
                }
            }
        });

        // xx
        pool->detach_task([this, &cx, q0, bins, ff_table, sinqd_table_xx] () {
            for (unsigned int q = q0; q < q0+bins; ++q) {
                auto xx = evaluate_xx_distance_profile(cx[q-q0]);
                double xx_sum = std::inner_product(xx.begin(), xx.end(), sinqd_table_xx->begin(q), 0.0);
                cache.intensity_profiles.xx[q-q0] += free_params.crho*free_params.crho*xx_sum*ff_table->index(form_factor::exv_bin, form_factor::exv_bin).evaluate(q);
            }
        });
    }

    if (cw_changed || cx_changed) {
        // wx
        pool->detach_task([this, &cx, q0, bins, ff_table, sinqd_table_ax] () {
            for (unsigned int q = q0; q < q0+bins; ++q) {
                auto wx = evaluate_wx_distance_profile(cx[q-q0]);
                double wx_sum = std::inner_product(wx.begin(), wx.end(), sinqd_table_ax->begin(q), 0.0);
                cache.intensity_profiles.wx[q-q0] += 2*free_params.crho*wx_sum*free_params.cw*ff_table->index(form_factor::water_bin, form_factor::exv_bin).evaluate(q);
            }
        });
    }
}
