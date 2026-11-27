// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <settings/SettingsValidation.h>

#include <settings/All.h>
#include <utility/Console.h>

#include <algorithm>
#include <cmath>
#include <string>

using namespace ausaxs;

void settings::validate_settings() {
    // check for exv fitting support: the simple models have no separate excluded volume to fit
    switch (settings::exv::exv_method.value) {
        case settings::exv::ExvMethod::None:
        case settings::exv::ExvMethod::Vacuum:
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

    {   // the grid-based excluded volume models place their dummies on every point_stride()'th grid cell, so a dummy
        // width that is not a whole number of cells is rounded. Only those models read it, so only they are warned
        switch (settings::exv::exv_method) {
            case settings::exv::ExvMethod::Grid:
            case settings::exv::ExvMethod::GridScalable:
            case settings::exv::ExvMethod::GridSurface:
            case settings::exv::ExvMethod::WAXSiS: {
                double ratio = settings::grid::exv::width/settings::grid::cell_width;
                if (1e-6 < std::abs(ratio - std::max(1., std::round(ratio)))) {
                    console::print_warning(
                        "Warning: The excluded volume width (" + std::to_string(settings::grid::exv::width) + " Å) is not a multiple of the grid cell width ("
                        + std::to_string(settings::grid::cell_width) + " Å). The dummy atoms will be spaced "
                        + std::to_string(std::max(1., std::round(ratio))*settings::grid::cell_width) + " Å apart instead."
                    );
                }
                break;
            }
            default: break;
        }
    }
}