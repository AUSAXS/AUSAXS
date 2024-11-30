/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <settings/SettingsValidation.h>
#include <settings/All.h>
#include <utility/Console.h>

using namespace ausaxs;

void settings::validate_settings() {
    switch (settings::hist::histogram_manager) {
        case settings::hist::HistogramManagerChoice::HistogramManagerMTFFAvg:
        case settings::hist::HistogramManagerChoice::HistogramManagerMTFFExplicit:
        case settings::hist::HistogramManagerChoice::FoXSManager:
        case settings::hist::HistogramManagerChoice::PepsiManager:
        case settings::hist::HistogramManagerChoice::CrysolManager:
        case settings::hist::HistogramManagerChoice::HistogramManagerMTFFGrid:
        case settings::hist::HistogramManagerChoice::HistogramManagerMTFFGridSurface:
            break;
        case settings::hist::HistogramManagerChoice::HistogramManagerMTFFGridScalableExv:
            if (settings::hist::histogram_manager == settings::hist::HistogramManagerChoice::HistogramManagerMTFFGridSurface) {
                // Use the faster basic Grid when the surface is not scaled anyway
                if (!settings::fit::fit_excluded_volume) {settings::hist::histogram_manager = settings::hist::HistogramManagerChoice::HistogramManagerMTFFGrid;}
            }
            break;
        case settings::hist::HistogramManagerChoice::PartialHistogramManager:
        case settings::hist::HistogramManagerChoice::PartialHistogramManagerMT:
            // check if a more efficient alternative is available
            if (settings::general::threads == 1) {
                console::print_warning("Warning: The chosen histogram manager is designed for multi-threading. Switching to single-threaded alternative.");
                if (settings::hist::histogram_manager == settings::hist::HistogramManagerChoice::PartialHistogramManager) {
                    settings::hist::histogram_manager = settings::hist::HistogramManagerChoice::HistogramManager;
                } else {
                    settings::hist::histogram_manager = settings::hist::HistogramManagerChoice::HistogramManagerMT;
                }
            }
            [[fallthrough]];
        default:
            // check for excluded volume fitting compatibility
            if (settings::fit::fit_excluded_volume) {
                console::print_warning("Warning: The chosen histogram manager does not support excluded volume fitting. Disabling excluded volume fitting.");
                settings::fit::fit_excluded_volume = false;
            }
    }

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
}