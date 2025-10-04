// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/intensity_calculator/CompositeDistanceHistogramFFGridScalableExv.h>
#include <settings/GridSettings.h>

using namespace ausaxs;
using namespace ausaxs::hist;

CompositeDistanceHistogramFFGridScalableExv::CompositeDistanceHistogramFFGridScalableExv(
    CompositeDistanceHistogramFFGrid&& res, 
    std::function<std::unique_ptr<CompositeDistanceHistogramFFGrid>(double)> eval_scaled_exv
) : CompositeDistanceHistogramFFGrid(std::move(res)), eval_scaled_exv(std::move(eval_scaled_exv)) {}

void CompositeDistanceHistogramFFGridScalableExv::apply_excluded_volume_scaling_factor(double k) {
    auto h = eval_scaled_exv(k);
    distance_profiles = std::move(h->distance_profiles);
    distance_axes = std::move(h->distance_axes);
    sinc_tables = std::move(h->sinc_tables);
    p = std::move(h->p);
    axis = std::move(h->axis);
    cache.sinqd.valid = false;
    auto V = std::pow(settings::grid::exv::width*k, 3);
    regenerate_ff_table(form_factor::ExvFormFactor(V));
}