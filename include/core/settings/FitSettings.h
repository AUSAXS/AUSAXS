// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <settings/ExportMacro.h>
#include <settings/SettingsHelper.h>

namespace ausaxs::settings {
    struct EXPORT fit {
        static Setting<bool> verbose;                 // Decides if the fitting process will be verbose.
        static Setting<int> N;                        // Number of points sampled when discretizing a model scattering curve
        static Setting<int> max_iterations;           // Maximum number of iterations in the fitting process
        static Setting<bool> fit_excluded_volume;     // Enable fitting of the excluded volume solvent.
        static Setting<bool> fit_solvent_density;     // Enable fitting of the solvent density for the excluded volume.
        static Setting<bool> fit_hydration;           // Enable fitting of the hydration shell.
        static Setting<bool> fit_atomic_debye_waller; // Enable fitting of the atomic form factor debye-waller factor.
        static Setting<bool> fit_exv_debye_waller;    // Enable fitting of the excluded volume form factor debye-waller factor.
    };
}