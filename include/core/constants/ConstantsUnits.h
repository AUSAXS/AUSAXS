// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <constants/ConstantsSI.h>
#include <math/ConstexprMath.h>

/**
 * @brief Relative units.
 * 
 * This namespace contains all the unit conversion constants used in this project. 
 * The basic units are for:
 *     mass: Dalton
 *     length: Å
 *     charge: e
 * 
 * Every factor here converts *out of* the internal units, so a published value enters them by
 * division: a density of 1.35 g/cm^3 is 1.35*unit::mL/unit::gm in Da/Å^3.
 */
namespace ausaxs::constants::unit {
    // Lengths are exact by definition of the SI prefixes, so they are written out directly.
    constexpr double cm = 1e-8; // Ångström --> cm
    constexpr double nm = 1e-1; // Ångström --> nm

    constexpr double mL = math::pow(unit::cm, 3); // Ångström^3 --> mL

    // The atomic mass unit is a *measured* constant, so it is derived from the single place it is
    // defined rather than restated in another base - the two spellings previously disagreed by
    // 6e-6 relative, which is larger than the constant's own uncertainty.
    constexpr double gm = SI::mass::u/SI::mass::gm; // Dalton --> grams
    constexpr double mg = SI::mass::u/SI::mass::mg; // Dalton --> mg
}
