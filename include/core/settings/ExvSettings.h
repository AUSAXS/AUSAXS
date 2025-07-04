// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <settings/ExportMacro.h>

namespace ausaxs::settings {
    struct EXPORT exv {
        enum class ExvMethod {
            // The simple method where the total charge of the excluded volume is subtracted from the atoms
            // An effective form factor of exp(-q*q) is used for all atoms
            Simple,

            // Equivalent to the simple method, except atoms now have their own individual form factors
            // The excluded volume contribution can be directly scaled by some factor
            Average,

            //? All following methods use individual atomic form factors

            // The excluded volume is modeled using a Gaussian spherical form factor at each atomic position,
            // with a radius dependent on the displaced volume of the atomic species.
            Fraser,

            // A regular grid is used to model a completely homogeneous excluded volume.
            // The grid cells share a single Gaussian spherical form factor, with a radius dependent on the cell volume.
            Grid,

            // Similar to the Grid method, except the cell width can be stretched to fit the excluded volume. 
            GridSurface,

            // Similar to the Grid method, except the occupied grid cells are marked as either interior or surface. 
            // The volume of only the surface cells can then be fitted in a similar manner to the Fraser method.
            GridScalable,

            // Mimic of the CRYSOL Fraser implementation.
            CRYSOL,

            // Mimic of the FoXS Fraser implementation.
            FoXS,

            // Mimic of the Pepsi Fraser implementation.
            Pepsi,

            // Mimic of the WAXSiS implementation, where the excluded volume covers the entire envelope.
            WAXSiS,
        };
        static ExvMethod exv_method; // The method used to model the excluded volume.
    };
}