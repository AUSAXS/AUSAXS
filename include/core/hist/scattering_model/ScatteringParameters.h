// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

namespace ausaxs::hist {
    /**
     * @brief The free parameters of the scattering models, at their default values. Each model only uses the ones it supports.
     */
    struct ScatteringParameters {
        double cw = 1;               // water density scaling factor
        double cx = 1;               // excluded volume scaling factor, method-dependent
        double crho = 1;             // solvent density scaling factor
        double DW_sigma_atomic = 0;  // atomic form factor debye-waller factor, zero for disabled
        double DW_sigma_exv = 0;     // excluded volume form factor debye-waller factor, zero for disabled
    };
}
