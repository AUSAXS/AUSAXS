// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <string>

namespace ausaxs::constants::fit {
    enum class Parameters {
        SCALING_WATER,
        SCALING_EXV,
        SCALING_RHO,
        DEBYE_WALLER_ATOMIC,
        DEBYE_WALLER_EXV
    };

    /**
     * @brief Get the short-form string representation of the given parameter.
     */
    inline std::string to_string(Parameters p) {
        switch (p) {
            case Parameters::SCALING_WATER:
                return "cw";
            case Parameters::SCALING_EXV:
                return "cx";
            case Parameters::SCALING_RHO:
                return "cr";
            case Parameters::DEBYE_WALLER_ATOMIC:
                return "Ba";
            case Parameters::DEBYE_WALLER_EXV:
                return "Bx";
        }
        return "";
    }
}