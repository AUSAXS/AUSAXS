// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <settings/ExportMacro.h>
#include <settings/SettingsHelper.h>

namespace ausaxs::settings {
    /// @brief Settings controlling the probe used in the scattering experiment.
    struct EXPORT scattering {
        /// @brief The available probes.
        enum class Radiation {
            XRay,   // Scattering from the electron clouds, described by the X-ray form factors.
            Neutron // Scattering from the nuclei, described by the coherent neutron scattering lengths.
        };

        // The probe used in the experiment. This decides which form factor tables are generated.
        // The effective charges of the atoms are derived from it when they are constructed, so it must be set before loading any structure.
        static detail::Setting<Radiation> radiation;
    };
}
