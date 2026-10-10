// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/intensity_calculator/CompositeDistanceHistogramFFGridBase.h>

#include <form_factor/ExvFormFactor.h>
#include <form_factor/lookup/FormFactorManager.h>
#include <settings/GridSettings.h>

using namespace ausaxs;
using namespace ausaxs::hist;
using namespace ausaxs::form_factor;

observer_ptr<const table::DebyeTable> CompositeDistanceHistogramFFGridBase::get_sinc_table_ax() const {
    return sinc_tables.ax.get_sinc_table();
}

observer_ptr<const table::DebyeTable> CompositeDistanceHistogramFFGridBase::get_sinc_table_xx() const {
    return sinc_tables.xx.get_sinc_table();
}

void CompositeDistanceHistogramFFGridBase::initialize_grid_axes(std::vector<double>&& d_axis_ax, std::vector<double>&& d_axis_xx) {
    distance_axes = {.xx=std::move(d_axis_xx), .ax=std::move(d_axis_ax)};
    sinc_tables.ax.set_d_axis(distance_axes.ax);
    sinc_tables.xx.set_d_axis(distance_axes.xx);
}

namespace {
    // Generate a form factor table for the grid-based calculations, using ffx as the excluded volume form factor.
    template<FormFactorType T>
    form_factor::lookup::table_t generate_ff_table(T&& ffx) {
        const auto* tables = form_factor::manager::get_active_product_tables();
        form_factor::manager::detail::profile_t ffx_profile;
        for (int q = 0; q < static_cast<int>(ffx_profile.size()); ++q) {ffx_profile[q] = ffx.evaluate(constants::axes::q_vals[q]);}

        // the atomic products are the same as those of the active tables; only the excluded volume row and column differ
        form_factor::lookup::table_t table = tables->raw_atomic_table;
        int n_active = form_factor::get_active_count();
        for (int i = 0; i < n_active; ++i) {
            auto atomic_profile = form_factor::manager::evaluate_amplitude(static_cast<form_factor_t>(tables->ff_indices[i]));
            table.index(i, form_factor::exv_bin) = FormFactorProduct(atomic_profile, ffx_profile);
            table.index(form_factor::exv_bin, i) = table.index(i, form_factor::exv_bin);
        }

        // must come last; the loop above overwrites this slot on its first iteration
        table.index(form_factor::exv_bin, form_factor::exv_bin) = FormFactorProduct(ffx_profile, ffx_profile);
        return table;
    }
}

template<FormFactorType T>
void CompositeDistanceHistogramFFGridBase::regenerate_ff_table(T&& ffx) {ff_table = generate_ff_table(std::forward<T>(ffx));}
template void CompositeDistanceHistogramFFGridBase::regenerate_ff_table(ExvFormFactor&&);

void CompositeDistanceHistogramFFGridBase::regenerate_ff_table() {
    regenerate_ff_table(ExvFormFactor(std::pow(settings::grid::exv::width, 3)));
}
