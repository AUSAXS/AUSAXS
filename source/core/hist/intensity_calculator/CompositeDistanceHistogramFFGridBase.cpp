// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/intensity_calculator/CompositeDistanceHistogramFFGridBase.h>
#include <form_factor/lookup/FormFactorManager.h>
#include <form_factor/lookup/NormalizedFormFactorProduct.h>
#include <form_factor/NormalizedFormFactor.h>
#include <form_factor/ExvFormFactor.h>
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

template<FormFactorType T>
void CompositeDistanceHistogramFFGridBase::regenerate_ff_table(T&& ffx) {ff_table = generate_ff_table(std::move(ffx));}
template void CompositeDistanceHistogramFFGridBase::regenerate_ff_table(ExvFormFactor&&);
template void CompositeDistanceHistogramFFGridBase::regenerate_ff_table(NormalizedFormFactor&&);

form_factor::lookup::table_t CompositeDistanceHistogramFFGridBase::generate_ff_table() {
    auto V = std::pow(settings::grid::exv::width, 3);
    return generate_ff_table(ExvFormFactor(V));
}

template<FormFactorType T>
form_factor::lookup::table_t CompositeDistanceHistogramFFGridBase::generate_ff_table(T&& ffx) {
    auto ff_indices = form_factor::manager::get_active_product_tables()->ff_indices;
    form_factor::lookup::table_t table;
    for (unsigned int i = 0; i < settings::form_factor::max_ff_types; ++i) {
        for (unsigned int j = 0; j < i; ++j) {
            table.index(i, j) = NormalizedFormFactorProduct(
                lookup::atomic::raw::get(static_cast<form_factor_t>(ff_indices[i])), 
                lookup::atomic::raw::get(static_cast<form_factor_t>(ff_indices[j]))
            );
            table.index(j, i) = table.index(i, j);
        }
        table.index(i, i) = NormalizedFormFactorProduct(
            lookup::atomic::raw::get(static_cast<form_factor_t>(ff_indices[i])), 
            lookup::atomic::raw::get(static_cast<form_factor_t>(ff_indices[i]))
        );

        table.index(i, form_factor::exv_bin) = NormalizedFormFactorProduct(
            lookup::atomic::raw::get(static_cast<form_factor_t>(ff_indices[i])), 
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
template form_factor::lookup::table_t CompositeDistanceHistogramFFGridBase::generate_ff_table(ExvFormFactor&&);
template form_factor::lookup::table_t CompositeDistanceHistogramFFGridBase::generate_ff_table(NormalizedFormFactor&&);
