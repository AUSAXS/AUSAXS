// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/intensity_calculator/CompositeDistanceHistogramFFGridSurface.h>
#include <form_factor/PrecalculatedFormFactorProduct.h>
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
) : hist::CompositeDistanceHistogramFFAvg(std::move(p_aa), std::move(p_aw), std::move(p_ww), std::move(p_tot_aa)),
    exv_distance_profiles{hist::Distribution1D(std::move(xx.interior)), hist::Distribution1D(std::move(xx.surface)), hist::Distribution1D(std::move(xx.cross)), 
                          hist::Distribution1D(std::move(wx.interior)), hist::Distribution1D(std::move(wx.surface)), hist::Distribution2D(std::move(ax.interior)), 
                          hist::Distribution2D(std::move(ax.surface))}
{
    initialize(p_tot_ax.get_weighted_axis(), p_tot_xx.get_weighted_axis());
}

Limit CompositeDistanceHistogramFFGridSurface::get_excluded_volume_scaling_factor_limits() const {
    return {0, 2};
}

const std::vector<double>& CompositeDistanceHistogramFFGridSurface::get_d_axis_xx() const {
    return distance_axes.xx;
}

const std::vector<double>& CompositeDistanceHistogramFFGridSurface::get_d_axis_ax() const {
    return distance_axes.ax;
}

const form_factor::storage::atomic::table_t& CompositeDistanceHistogramFFGridSurface::get_ff_table() const {
    return CompositeDistanceHistogramFFGrid::ff_table;
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

observer_ptr<const table::DebyeTable> CompositeDistanceHistogramFFGridSurface::get_sinc_table_ax() const {
    return sinc_tables.ax.get_sinc_table();
}

observer_ptr<const table::DebyeTable> CompositeDistanceHistogramFFGridSurface::get_sinc_table_xx() const {
    return sinc_tables.xx.get_sinc_table();
}

void CompositeDistanceHistogramFFGridSurface::initialize(std::vector<double>&& d_axis_ax, std::vector<double>&& d_axis_xx) {
    static bool initialized = false;
    if (!initialized) {
        CompositeDistanceHistogramFFGrid::ff_table = CompositeDistanceHistogramFFGrid::generate_ff_table();
        initialized = true;
    }

    this->distance_axes = {.xx=std::move(d_axis_xx), .ax=std::move(d_axis_ax)};
    sinc_tables.ax.set_d_axis(this->distance_axes.ax);
    sinc_tables.xx.set_d_axis(this->distance_axes.xx);

    // fix the aa counts to also contain the exv contributions
    auto xx = evaluate_xx_distance_profile(1);
    auto wx = evaluate_wx_distance_profile(1);
    auto ax = evaluate_ax_distance_profile(1);

    auto& aa = CompositeDistanceHistogramFFAvgBase::get_aa_counts_by_ff();
    for (unsigned int ff = 0; ff < form_factor::get_count_without_excluded_volume(); ++ff) {
        std::transform(aa.begin(ff, form_factor::exv_bin), aa.end(ff, form_factor::exv_bin), ax.begin(ff), aa.begin(ff, form_factor::exv_bin), std::plus<double>());
    }
    std::transform(aa.begin(form_factor::exv_bin, form_factor::exv_bin), aa.end(form_factor::exv_bin, form_factor::exv_bin), xx.begin(), aa.begin(form_factor::exv_bin, form_factor::exv_bin), std::plus<double>());

    auto& aw = CompositeDistanceHistogramFFAvgBase::get_aw_counts_by_ff();
    std::transform(aw.begin(form_factor::exv_bin), aw.end(form_factor::exv_bin), wx.begin(), aw.begin(form_factor::exv_bin), std::plus<double>());
}

void CompositeDistanceHistogramFFGridSurface::cache_refresh_sinqd() const {
    auto pool = utility::multi_threading::get_global_pool();
    const auto& sinqd_table = sinc_table.get_sinc_table();

    Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
    unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin);

    if (cache.sinqd.aa.empty()) {
        cache.sinqd.aa = container::Container3D<double>(form_factor::get_count(), form_factor::get_count(), debye_axis.bins);
        cache.sinqd.aw = container::Container2D<double>(form_factor::get_count(), debye_axis.bins);
        cache.sinqd.ww = container::Container1D<double>(debye_axis.bins);
    }

    for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
        for (unsigned int ff2 = 0; ff2 < form_factor::get_count_without_excluded_volume(); ++ff2) {
            pool->detach_task([this, q0, bins=debye_axis.bins, ff1, ff2, sinqd_table] () {
                for (unsigned int q = q0; q < q0+bins; ++q) {
                    cache.sinqd.aa.index(ff1, ff2, q-q0) = std::inner_product(distance_profiles.aa.begin(ff1, ff2), distance_profiles.aa.end(ff1, ff2), sinqd_table->begin(q), 0.0);
                }
            });
        }
        pool->detach_task([this, q0, bins=debye_axis.bins, ff1, sinqd_table] () {
            for (unsigned int q = q0; q < q0+bins; ++q) {
                cache.sinqd.aw.index(ff1, q-q0) = std::inner_product(distance_profiles.aw.begin(ff1), distance_profiles.aw.end(ff1), sinqd_table->begin(q), 0.0);
            }
        });
    }
    pool->detach_task([&] () {
        for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
            cache.sinqd.ww.index(q-q0) = std::inner_product(distance_profiles.ww.begin(), distance_profiles.ww.end(), sinqd_table->begin(q), 0.0);
        }
    });
    cache.sinqd.valid = true;
    pool->wait();
}

void CompositeDistanceHistogramFFGridSurface::cache_refresh_intensity_profiles(bool sinqd_changed, bool cw_changed, bool cx_changed) const {
    auto pool = utility::multi_threading::get_global_pool();
    const auto& ff_table = get_ff_table();
    auto sinqd_table_ax = get_sinc_table_ax();
    auto sinqd_table_xx = get_sinc_table_xx();

    Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
    unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

    if (sinqd_changed) {
        this->cache.intensity_profiles.aa = std::vector<double>(debye_axis.bins, 0);
    }
    if (cw_changed) {
        this->cache.intensity_profiles.aw = std::vector<double>(debye_axis.bins, 0);
        this->cache.intensity_profiles.ww = std::vector<double>(debye_axis.bins, 0);
    }
    if (cx_changed) {
        this->cache.intensity_profiles.ax = std::vector<double>(debye_axis.bins, 0);
        this->cache.intensity_profiles.xx = std::vector<double>(debye_axis.bins, 0);
    }
    if (cw_changed || cx_changed) {
        this->cache.intensity_profiles.wx = std::vector<double>(debye_axis.bins, 0);
    }

    std::vector<double> cx(debye_axis.bins, 0);
    for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {cx[q-q0] = exv_factor(constants::axes::q_vals[q]);}

    if (sinqd_changed) {
        // aa
        pool->detach_task([&] () {
            for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
                for (unsigned int ff2 = 0; ff2 < form_factor::get_count_without_excluded_volume(); ++ff2) {
                    for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
                        this->cache.intensity_profiles.aa[q-q0] += 
                            this->cache.sinqd.aa.index(ff1, ff2, q-q0)*ff_table.index(ff1, ff2).evaluate(q);
                    }
                }
            }
        });
    }

    if (cx_changed) {
        // ax
        pool->detach_task([&] () {
            for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
                auto ax = evaluate_ax_distance_profile(cx[q-q0]);
                for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
                    double ax_sum = std::inner_product(ax.begin(ff1), ax.end(ff1), sinqd_table_ax->begin(q), 0.0);
                    this->cache.intensity_profiles.ax[q-q0] += 
                        2*this->free_params.crho*ax_sum*ff_table.index(ff1, form_factor::exv_bin).evaluate(q);
                }
            }
        });

        // xx
        pool->detach_task([&] () {
            for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
                auto xx = evaluate_xx_distance_profile(cx[q-q0]);
                double xx_sum = std::inner_product(xx.begin(), xx.end(), sinqd_table_xx->begin(q), 0.0);
                this->cache.intensity_profiles.xx[q-q0] += 
                    this->free_params.crho*this->free_params.crho*xx_sum*ff_table.index(form_factor::exv_bin, form_factor::exv_bin).evaluate(q);
            }
        });
    }

    if (cw_changed) {
        // aw
        pool->detach_task([&] () {
            for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
                for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
                    this->cache.intensity_profiles.aw[q-q0] += 
                        2*this->free_params.cw*this->cache.sinqd.aw.index(ff1, q-q0)
                        *ff_table.index(ff1, form_factor::water_bin).evaluate(q);
                }
            }
        });

        // ww
        pool->detach_task([&] () {
            for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
                this->cache.intensity_profiles.ww[q-q0] += 
                    this->free_params.cw*this->free_params.cw*this->cache.sinqd.ww.index(q-q0)
                    *ff_table.index(form_factor::water_bin, form_factor::water_bin).evaluate(q);
            }
        });
    }

    if (cw_changed || cx_changed) {
        // wx
        pool->detach_task([&] () {
            for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
                auto wx = evaluate_wx_distance_profile(cx[q-q0]);
                double wx_sum = std::inner_product(wx.begin(), wx.end(), sinqd_table_ax->begin(q), 0.0);
                this->cache.intensity_profiles.wx[q-q0] += 
                    2*this->free_params.crho*wx_sum*this->free_params.cw*ff_table.index(form_factor::water_bin, form_factor::exv_bin).evaluate(q);
            }
        });
    }
    this->cache.intensity_profiles.cached_cx = this->free_params.cx;
    this->cache.intensity_profiles.cached_cw = this->free_params.cw;
    this->cache.intensity_profiles.cached_crho = this->free_params.crho;
    pool->wait();
}