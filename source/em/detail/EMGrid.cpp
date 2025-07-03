// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <em/detail/EMGrid.h>
#include <settings/GridSettings.h>

using namespace ausaxs::em::grid;

double EMGrid::get_atomic_radius(form_factor::form_factor_t) const {
    return settings::grid::min_exv_radius;
}