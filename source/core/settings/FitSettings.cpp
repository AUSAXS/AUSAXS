// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

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
    settings::io::SettingSection fit_section("Fit", {
        settings::io::create(fit::verbose, "fit-verbose"),
        settings::io::create(fit::N, "N"),
        settings::io::create(fit::max_iterations, "max_iterations"),
        settings::io::create(fit::fit_excluded_volume, "excluded_volume"),
        settings::io::create(fit::fit_solvent_density, "solvent_density"),
        settings::io::create(fit::fit_hydration, "hydration")
    });
}