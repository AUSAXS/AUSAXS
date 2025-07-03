// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <settings/SettingsValidation.h>
#include <settings/All.h>
#include <constants/ConstantsAxes.h>
#include <utility/Console.h>
#include <utility/Logging.h>

using namespace ausaxs;

void settings::validate_settings() {
    // check for exv fitting support
    switch (settings::hist::get_histogram_manager()) {
        // the following managers do not support exv fitting
        case settings::hist::HistogramManagerChoice::HistogramManager:
        case settings::hist::HistogramManagerChoice::HistogramManagerMT:
        case settings::hist::HistogramManagerChoice::HistogramSymmetryManagerMT:
        case settings::hist::HistogramManagerChoice::PartialHistogramManager:
        case settings::hist::HistogramManagerChoice::PartialHistogramManagerMT:
        case settings::hist::HistogramManagerChoice::PartialHistogramSymmetryManagerMT:
            if (settings::fit::fit_excluded_volume) {
                console::print_warning("Warning: The chosen histogram manager does not support excluded volume fitting. Disabling excluded volume fitting.");
                settings::fit::fit_excluded_volume = false;
            }

        // we explicitly write each case to ensure we will get a compiler warning for new managers in the future
        case settings::hist::HistogramManagerChoice::HistogramManagerMTFFAvg:
        case settings::hist::HistogramManagerChoice::HistogramManagerMTFFExplicit:
        case settings::hist::HistogramManagerChoice::FoXSManager:
        case settings::hist::HistogramManagerChoice::PepsiManager:
        case settings::hist::HistogramManagerChoice::CrysolManager:
        case settings::hist::HistogramManagerChoice::HistogramManagerMTFFGrid:
        case settings::hist::HistogramManagerChoice::HistogramManagerMTFFGridScalableExv:
        case settings::hist::HistogramManagerChoice::HistogramManagerMTFFGridSurface:
        case settings::hist::HistogramManagerChoice::Count:
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
        if (std::abs(grid_ratio - int(grid_ratio))) {
            console::print_warning("Warning: The grid cell width is not a multiple of the excluded volume radius. This may lead to artifacts in the excluded volume calculation.");
        }
    }

    {   // if the user wants to keep hydrogens, they should be treated explicitly
        if (settings::general::keep_hydrogens) {
            settings::molecule::implicit_hydrogens = false;
        } else {
            settings::molecule::implicit_hydrogens = true;
        }
    }

    {   // check for q-range compatibility
        if (settings::axes::qmin < constants::axes::q_axis.min) {
            console::print_warning("Warning: qmin is smaller than the minimum q value. Setting qmin to the minimum q value (" + std::to_string(constants::axes::q_axis.min) + ").");
            settings::axes::qmin = constants::axes::q_axis.min;
        }

        if (constants::axes::q_axis.max < settings::axes::qmax) {
            console::print_warning("Warning: qmax is larger than the maximum q value. Setting qmax to the maximum q value (" + std::to_string(constants::axes::q_axis.max) + ").");
            settings::axes::qmax = constants::axes::q_axis.max;
        }

        if (settings::axes::qmax < settings::axes::qmin) {
            console::print_warning("Warning: qmax is smaller than qmin. Resetting to default values (0, 0.5).");
            settings::axes::qmin = 0;
            settings::axes::qmax = 0.5;
        }
    }
}