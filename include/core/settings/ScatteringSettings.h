// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <settings/ExportMacro.h>
#include <settings/SettingsHelper.h>

namespace ausaxs::settings {
    struct EXPORT scattering {
        enum class Radiation {
            XRay,
            Neutron
        };

        static Setting<Radiation> radiation;
    };
}
