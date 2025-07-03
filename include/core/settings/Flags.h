// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

// Internal boolean flags for signalling state. 
// This is primarily used to communicate between different parts of the program.
namespace ausaxs::settings {
    struct flags {
        static bool data_rebin;             // Whether the data has been rebinned.
        static char last_parsed_unit;       // The unit parsed from the latest file. 
        static bool init_histogram_manager; // Whether to initialize the histogram manager when a Molecule is created.
    };
}