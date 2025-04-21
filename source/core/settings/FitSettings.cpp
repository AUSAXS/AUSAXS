/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <settings/FitSettings.h>
#include <settings/SettingsIORegistry.h>

using namespace ausaxs;

bool settings::fit::verbose = false;
unsigned int settings::fit::N = 100;
unsigned int settings::fit::max_iterations = 100;
bool settings::fit::fit_excluded_volume = false;
bool settings::fit::fit_solvent_density = false;
bool settings::fit::fit_hydration = true;
bool settings::fit::fit_atomic_debye_waller = false;
bool settings::fit::fit_exv_debye_waller = false;

namespace ausaxs::settings::io {
    settings::io::SettingSection fit_section("General", {
        settings::io::create(fit::verbose, "fit-verbose"),
        settings::io::create(fit::N, "N"),
        settings::io::create(fit::max_iterations, "max_iterations"),
        settings::io::create(fit::fit_excluded_volume, "fit_excluded_volume"),
        settings::io::create(fit::fit_solvent_density, "fit_solvent_density"),
        settings::io::create(fit::fit_hydration, "fit_hydration")
    });
}