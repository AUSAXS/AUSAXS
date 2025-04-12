#pragma once

#include <string>

namespace ausaxs::settings::md {
    extern std::string gmx_path;                    // The path to the GROMACS installation.
    extern std::string gmx_source_path;             // The path to the GROMACS source installation.
    extern std::string buffer_path;                 // The path to the buffer directory.

    std::string gmx_top_path();
}