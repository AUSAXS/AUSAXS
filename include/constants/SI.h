#pragma once

/**
 * @brief Absolute units.
 * 
 * This namespace contains all the absolute unit conversion constants. 
 */
namespace constants::SI {
    namespace mass {
        constexpr double kg = 1;
        constexpr double gm = 1e-3;
        constexpr double mg = 1e-6;
        constexpr double u = 1.66053*1e-27;
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