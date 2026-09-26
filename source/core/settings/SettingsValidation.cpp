// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <settings/SettingsValidation.h>

#include <hist/detail/SimpleExvModel.h>
#include <settings/All.h>
#include <utility/Console.h>

using namespace ausaxs;

void settings::validate_settings() {
    // check for exv fitting support: the simple models have no separate excluded volume to fit
    switch (settings::exv::exv_method) {
        case settings::exv::ExvMethod::None:
            ausaxs::hist::detail::SimpleExvModel::disable(); // no excluded volume at all, so not even the effective charges of the simple model
            [[fallthrough]];
        case settings::exv::ExvMethod::Simple:
            if (settings::fit::fit_excluded_volume) {
                console::print_warning("Warning: The chosen excluded volume model does not support excluded volume fitting. Disabling excluded volume fitting.");
                settings::fit::fit_excluded_volume = false;
            }
            break;

        // we explicitly write each case to ensure we will get a compiler warning for new models in the future
        case settings::exv::ExvMethod::Average:
        case settings::exv::ExvMethod::Fraser:
        case settings::exv::ExvMethod::Grid:
        case settings::exv::ExvMethod::GridSurface:
        case settings::exv::ExvMethod::GridScalable:
        case settings::exv::ExvMethod::CRYSOL:
        case settings::exv::ExvMethod::FoXS:
        case settings::exv::ExvMethod::Pepsi:
        case settings::exv::ExvMethod::WAXSiS:
            break;
    }

    // if the pepsi mimic exv method is used, also match the cell widths and hydration strategy to theirs
    switch (settings::hydrate::hydration_strategy) {
        case settings::hydrate::HydrationStrategy::PepsiStrategy:
            if (settings::grid::cell_width < 3) {
                console::print_warning("Warning: The Pepsi hydration method requires a specific set of grid options. Setting grid width to 3Å and all atomic radii to 3Å.");
                settings::grid::cell_width = 3;
                settings::grid::min_exv_radius = 3;
            }
            break;
        default:
            break;
    }

    {   // check for grid cell width compatibility
        double grid_ratio = settings::grid::exv::width/settings::grid::cell_width;
        if (std::abs(grid_ratio - int(grid_ratio)) != 0.0) {
            console::print_warning("Warning: The grid cell width is not a multiple of the excluded volume radius. This may lead to artifacts in the excluded volume calculation.");
        }
    }

    {   // if the user wants to keep hydrogens, they should be treated explicitly
        settings::molecule::implicit_hydrogens = !settings::general::keep_hydrogens;
    }
}