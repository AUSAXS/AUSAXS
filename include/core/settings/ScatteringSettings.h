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

        // The probe. It selects which of the xray and neutron settings below are in effect.
        static detail::Setting<Radiation> radiation;

        // Settings which only apply to X-rays.
        struct xray {
            // The electron density of the bulk solvent in e/Å^3, which scales every excluded volume. The default is pure H2O.
            static detail::Setting<double> solvent_density;
        };

        // Settings which only apply to neutrons.
        struct neutron {
            // The coherent scattering length density of the bulk solvent in fm/Å^3, which scales every excluded volume. The default is pure D2O.
            static detail::Setting<double> solvent_density;

            // The fraction of the exchangeable implicit hydrogens (those bound to N, O or S) that has exchanged with the solvent. The default is full exchange.
            static detail::Setting<double> exchange_efficiency;

            // Set solvent_density to the density of an H2O/D2O mixture with the given D2O fraction, between 0 and 1.
            static void set_d2o_fraction(double d2o_fraction);
        };
    };
}
