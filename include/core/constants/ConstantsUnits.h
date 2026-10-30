// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <constants/ConstantsSI.h>
#include <math/ConstexprMath.h>

/**
 * @brief Relative units.
 * 
 * The internal unit system is:
 *     mass: Dalton
 *     length: Å
 *     charge: e
 * 
 * Multiplying by a unit::internal factor brings a value into it, and multiplying by the matching
 * unit::external factor takes it back out.
 */
namespace ausaxs::constants::unit {
    namespace internal {
        constexpr double cm = 1e8;
        constexpr double nm = 1e1;
        constexpr double mL = math::pow(cm, 3);
        constexpr double gm = SI::mass::gm/SI::mass::u;
        constexpr double mg = SI::mass::mg/SI::mass::u;
    }

    namespace external {
        constexpr double cm = 1/internal::cm;
        constexpr double nm = 1/internal::nm;
        constexpr double mL = 1/internal::mL;
        constexpr double gm = 1/internal::gm;
        constexpr double mg = 1/internal::mg;
    }
}
