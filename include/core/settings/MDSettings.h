#pragma once

#include <settings/ExportMacro.h>

#include <string>

namespace ausaxs::settings {
    struct EXPORT md {
        static std::string gmx_path;                    // The path to the GROMACS installation.
        static std::string gmx_source_path;             // The path to the GROMACS source installation.
        static std::string buffer_path;                 // The path to the buffer directory.
    
        std::string gmx_top_path();
    };
}