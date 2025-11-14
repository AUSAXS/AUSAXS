// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

namespace ausaxs::settings {
    // Internal data variables for signalling state. 
    // Do not manually modify these in your scripts or bad things will happen. 
    struct flags {
        static bool data_rebin;             // Whether the data has been rebinned.
        static char last_parsed_unit;       // The unit parsed from the latest file. 
        static bool init_histogram_manager; // Whether to initialize the histogram manager when a Molecule is created.
        static bool custom_bin_width;       // Whether a custom bin width is being used for the distance histogram.
        static double inv_bin_width;        // The inverse of the bin width for the distance histogram.
    };
}