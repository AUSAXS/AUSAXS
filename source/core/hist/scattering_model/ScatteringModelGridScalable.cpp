// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <form_factor/ExvFormFactor.h>
#include <hist/scattering_model/ScatteringModelGridScalable.h>

#include <settings/GridSettings.h>

using namespace ausaxs;
using namespace ausaxs::hist;

ScatteringModelGridScalable::ScatteringModelGridScalable(
    ScatteringModelGrid&& res, 
    std::function<std::unique_ptr<ScatteringModelGrid>(double)> eval_scaled_exv
) : ScatteringModelGrid(std::move(res)), eval_scaled_exv(std::move(eval_scaled_exv)) {}

void ScatteringModelGridScalable::apply_excluded_volume_scaling_factor(double k) {
    auto h = eval_scaled_exv(k);
    distance_profiles = std::move(h->distance_profiles);
    distance_axes = std::move(h->distance_axes);
    sinc_tables = std::move(h->sinc_tables);
    p = std::move(h->p);
    axis = h->axis;
    cache.sinqd.valid = false;
    auto V = std::pow(settings::grid::exv::width*k, 3);
    regenerate_ff_table(form_factor::factory::construct_exv_form_factor(V));
}