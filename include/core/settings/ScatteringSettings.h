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

        static detail::Setting<Radiation> radiation;

        // The electron density of the bulk solvent in e/Å^3, which scales every excluded volume for X-rays. The default is pure H2O.
        static detail::Setting<double> xray_solvent_density;

        // The coherent scattering length density of the bulk solvent in fm/Å^3, which scales every excluded volume for neutrons. The default is pure H2O.
        // Use form_factor::neutron::solvent_density for an H2O/D2O mixture.
        static detail::Setting<double> neutron_solvent_density;
    };
}
