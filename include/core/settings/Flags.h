#pragma once

// Namespace for internal boolean flags to signal state. 
// This is primarily used to communicate between different parts of the program.
namespace ausaxs::settings {
    struct flags {
        static bool data_rebin;         // Whether the data has been rebinned.
        static char last_parsed_unit;   // The unit parsed from the latest file. 
    };
}