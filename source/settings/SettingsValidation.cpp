/*
This software is distributed under the GNU General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <settings/SettingsValidation.h>
#include <settings/All.h>
#include <utility/Console.h>

void settings::validate_settings() {
    switch(settings::hist::histogram_manager) {
        case settings::hist::HistogramManagerChoice::HistogramManagerMTFFAvg:
        case settings::hist::HistogramManagerChoice::HistogramManagerMTFFExplicit:
        case settings::hist::HistogramManagerChoice::HistogramManagerMTFFGrid:
            // check for effective charge compatibility
            if (settings::molecule::use_effective_charge == true) {
                console::print_warning("Warning: The chosen histogram manager does not support using an effective charge approximation.");
                settings::molecule::use_effective_charge = false;
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
            if (settings::hist::fit_excluded_volume) {
                console::print_warning("Warning: The chosen histogram manager does not support excluded volume fitting.");
                settings::hist::fit_excluded_volume = false;
            }
    }
}