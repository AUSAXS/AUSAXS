// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

/**
 * @brief Absolute units.
 * 
 * This namespace contains all the absolute unit conversion constants. 
 * Each constant is the factor that takes a value *to* SI, not a value that is already in SI:
 * SI::length::A is "one Ångström, expressed in metres", so a length in Å becomes metres by
 * multiplication. See constants::unit for the conversions out of the internal unit system.
 */
namespace ausaxs::constants::SI {
    namespace mass {
        constexpr double kg = 1;
        constexpr double gm = 1e-3;
        constexpr double mg = 1e-6;
        constexpr double u = 1.66053906892e-27; // CODATA 2022: 1.66053906892(52)e-27 kg
    }

    namespace length {
        constexpr double m = 1;
        constexpr double cm = 1e-2;
        constexpr double nm = 1e-9;
        constexpr double A = 1e-10; // Ångström
    }

    namespace volume {
        constexpr double A3 = 1e-30; // Ångström^3
        constexpr double cm3 = 1e-6;
    }
}
