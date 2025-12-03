// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/intensity_calculator/CompositeDistanceHistogramFFGrid.h>
#include <form_factor/lookup/NormalizedFormFactorProduct.h>
#include <form_factor/NormalizedFormFactor.h>
#include <form_factor/ExvFormFactor.h>
#include <table/ArrayDebyeTable.h>
#include <settings/GridSettings.h>
#include <settings/HistogramSettings.h>
#include <dataset/SimpleDataset.h>
#include <utility/Exceptions.h>
#include <utility/MultiThreading.h>

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
) : hist::CompositeDistanceHistogramFFAvg(std::move(p_aa), std::move(p_aw), std::move(p_ww), std::move(p_tot_aa)) {
    initialize(p_tot_ax.get_weighted_axis(), p_tot_xx.get_weighted_axis());
}

template<FormFactorType T>
void CompositeDistanceHistogramFFGrid::regenerate_ff_table(T&& ffx) {ff_table = generate_ff_table(std::move(ffx));}
template void CompositeDistanceHistogramFFGrid::regenerate_ff_table(ExvFormFactor&&);
template void CompositeDistanceHistogramFFGrid::regenerate_ff_table(NormalizedFormFactor&&);

double CompositeDistanceHistogramFFGrid::exv_factor(double) const {
    return free_params.cx;
}

form_factor::lookup::atomic::table_t CompositeDistanceHistogramFFGrid::generate_ff_table() {
    auto V = std::pow(settings::grid::exv::width, 3);
    return generate_ff_table(ExvFormFactor(V));
}

template<FormFactorType T>
form_factor::lookup::atomic::table_t CompositeDistanceHistogramFFGrid::generate_ff_table(T&& ffx) {
    form_factor::lookup::atomic::table_t table;
    for (unsigned int i = 0; i < form_factor::get_count_without_excluded_volume(); ++i) {
        for (unsigned int j = 0; j < i; ++j) {
            table.index(i, j) = NormalizedFormFactorProduct(
                lookup::atomic::raw::get(static_cast<form_factor_t>(i)), 
                lookup::atomic::raw::get(static_cast<form_factor_t>(j))
            );
            table.index(j, i) = table.index(i, j);
        }
        table.index(i, i) = NormalizedFormFactorProduct(
            lookup::atomic::raw::get(static_cast<form_factor_t>(i)), 
            lookup::atomic::raw::get(static_cast<form_factor_t>(i))
        );

        table.index(i, form_factor::exv_bin) = NormalizedFormFactorProduct(
            lookup::atomic::raw::get(static_cast<form_factor_t>(i)), 
            ffx
        );
        table.index(form_factor::exv_bin, i) = table.index(i, form_factor::exv_bin);
        table.index(form_factor::exv_bin, form_factor::exv_bin) = NormalizedFormFactorProduct(
            ffx, 
            ffx
        );
    }
    return table;
}
template form_factor::lookup::atomic::table_t CompositeDistanceHistogramFFGrid::generate_ff_table(ExvFormFactor&&);
template form_factor::lookup::atomic::table_t CompositeDistanceHistogramFFGrid::generate_ff_table(NormalizedFormFactor&&);

void CompositeDistanceHistogramFFGrid::cache_refresh_sinqd() const {
    auto pool = utility::multi_threading::get_global_pool();
    const auto& sinqd_table_aa = sinc_table.get_sinc_table();
    const auto& sinqd_table_ax = get_sinc_table_ax();
    const auto& sinqd_table_xx = get_sinc_table_xx();

    Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
    unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin);

    if (cache.sinqd.aa.empty()) {
        cache.sinqd.aa = container::Container3D<double>(form_factor::get_count(), form_factor::get_count(), debye_axis.bins);
        cache.sinqd.ax = container::Container2D<double>(form_factor::get_count(), debye_axis.bins);
        cache.sinqd.aw = container::Container2D<double>(form_factor::get_count(), debye_axis.bins);
        cache.sinqd.xx = container::Container1D<double>(debye_axis.bins);
        cache.sinqd.wx = container::Container1D<double>(debye_axis.bins);
        cache.sinqd.ww = container::Container1D<double>(debye_axis.bins);
    }

    for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
        for (unsigned int ff2 = 0; ff2 < form_factor::get_count_without_excluded_volume(); ++ff2) {
            pool->detach_task([this, q0, bins=debye_axis.bins, ff1, ff2, sinqd_table_aa] () {
                for (unsigned int q = q0; q < q0+bins; ++q) {
                    cache.sinqd.aa.index(ff1, ff2, q-q0) = std::inner_product(distance_profiles.aa.begin(ff1, ff2), distance_profiles.aa.end(ff1, ff2), sinqd_table_aa->begin(q), 0.0);
                }
            });
        }
        pool->detach_task([this, q0, bins=debye_axis.bins, ff1, sinqd_table_aa, sinqd_table_ax] () {
            for (unsigned int q = q0; q < q0+bins; ++q) {
                cache.sinqd.ax.index(ff1, q-q0) = std::inner_product(distance_profiles.aa.begin(ff1, form_factor::exv_bin), distance_profiles.aa.end(ff1, form_factor::exv_bin), sinqd_table_ax->begin(q), 0.0);
                cache.sinqd.aw.index(ff1, q-q0) = std::inner_product(distance_profiles.aw.begin(ff1), distance_profiles.aw.end(ff1), sinqd_table_aa->begin(q), 0.0);
            }
        });
    }
    pool->detach_task([&] () {
        for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
            cache.sinqd.xx.index(q-q0) = std::inner_product(distance_profiles.aa.begin(form_factor::exv_bin, form_factor::exv_bin), distance_profiles.aa.end(form_factor::exv_bin, form_factor::exv_bin), sinqd_table_xx->begin(q), 0.0);
            cache.sinqd.wx.index(q-q0) = std::inner_product(distance_profiles.aw.begin(form_factor::exv_bin), distance_profiles.aw.end(form_factor::exv_bin), sinqd_table_ax->begin(q), 0.0);
            cache.sinqd.ww.index(q-q0) = std::inner_product(distance_profiles.ww.begin(), distance_profiles.ww.end(), sinqd_table_aa->begin(q), 0.0);
        }
    });
    cache.sinqd.valid = true;
    pool->wait();
}

observer_ptr<const table::DebyeTable> CompositeDistanceHistogramFFGrid::get_sinc_table_ax() const {
    return sinc_tables.ax.get_sinc_table();
}

observer_ptr<const table::DebyeTable> CompositeDistanceHistogramFFGrid::get_sinc_table_xx() const {
    return sinc_tables.xx.get_sinc_table();
}

void CompositeDistanceHistogramFFGrid::initialize(std::vector<double>&& d_axis_ax, std::vector<double>&& d_axis_xx) {
    static bool init_table = false;
    if (!init_table) {
        ff_table = generate_ff_table();
        init_table = true;
    }

    this->distance_axes = {.xx=std::move(d_axis_xx), .ax=std::move(d_axis_ax)};
    sinc_tables.ax.set_d_axis(this->distance_axes.ax);
    sinc_tables.xx.set_d_axis(this->distance_axes.xx);
}