// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <settings/ExportMacro.h>

namespace ausaxs::settings {
    // Internal data variables for signalling state. 
    // Do not manually modify these in your scripts or bad things will happen. 
    struct EXPORT internal_state {
        static bool custom_bin_width;       // Whether a custom bin width is being used for the distance histogram.
        static double inv_bin_width;        // The inverse of the bin width for the distance histogram.
        static bool prefer_partial_manager; // Whether to prefer a partial histogram manager if one is available for the chosen excluded volume method.
        static bool allow_decorrelate_atom_order; // Whether the distance calculators may permute the atom order before accumulating a histogram.
    };
}
