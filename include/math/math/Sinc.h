// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <math/ConstexprMath.h>

#include <cmath>
#include <type_traits>

namespace ausaxs::math {
    /**
     * @brief The unnormalized sinc function, sin(x)/x, with sinc(0) = 1.
     *        Usable in constant expressions.
     */
    constexpr double sinc(double x) {
        // sin(x)/x stays accurate for small x, so the series only removes the 0/0 at the origin.
        // Its first omitted term, x^6/5040, is below double precision for |x| < 1e-3.
        if (-1e-3 < x && x < 1e-3) {
            double x2 = x*x;
            return 1 - x2/6 + x2*x2/120;
        }
        if (std::is_constant_evaluated()) {return constexpr_math::sin(x)/x;}
        return std::sin(x)/x;
    }
}
