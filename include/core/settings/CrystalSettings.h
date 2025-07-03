// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <settings/ExportMacro.h>

namespace ausaxs::settings {
    struct EXPORT crystal {
        static unsigned int h;        // The maximum Miller index along the x direction.
        static unsigned int k;        // The maximum Miller index along the y direction.
        static unsigned int l;        // The maximum Miller index along the z direction.
    
        static double max_q;          // The maximum length of the Miller indices. 
        static double grid_expansion; // The factor by which the grid is expanded when loading a pdb structure.

        struct reduced {
            static double basis_q;    // The maximum q value for which the basis is generated.
        };

        struct detail {
            static bool use_checkpointing; // Whether to use checkpointing during the calculation. 
        };

        // The choice of Miller index generation algorithm. See the documentation for the individual classes for more information.
        enum class MillerGenerationChoice {
            All, // Generates all miller indices within the specified range. 
            Reduced, // Generates a subset of the miller indices within the specified range using a reduced basis.
            Fibonacci // Similar to Reduced, but uses a Fibonacci sphere to make the basis more uniform. Experimental and not recommended for general usage.
        };
        static MillerGenerationChoice miller_generation_strategy;
    };
}