// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <gcem.hpp>

#include <cmath>
#include <type_traits>
#include <version>

// Math functions usable in constant expressions. Constant evaluation uses the standard library where
// it provides constexpr math functions, otherwise the third-party gcem library as a fallback.
// At runtime they always call the standard library, since gcem is not meant for runtime use.
namespace ausaxs::math {
    namespace detail {
        #if defined __cpp_lib_constexpr_cmath
            namespace backend = std;
        #else
            namespace backend = ::gcem;
        #endif
    }

    constexpr double sin(double x) {
        if (std::is_constant_evaluated()) {return detail::backend::sin(x);}
        return std::sin(x);
    }

    constexpr double sinc(double x) {
        if (-1e-3 < x && x < 1e-3) {
            double x2 = x*x;
            return 1 - x2/6 + x2*x2/120;
        }
        return sin(x)/x;
    }

    constexpr double exp(double x) {
        if (std::is_constant_evaluated()) {return detail::backend::exp(x);}
        return std::exp(x);
    }

    // the exponent type is kept, so gcem can raise integral exponents exactly by repeated squaring
    template<typename T> requires std::is_arithmetic_v<T>
    constexpr double pow(double base, T exponent) {
        if (std::is_constant_evaluated()) {return detail::backend::pow(base, exponent);}
        return std::pow(base, exponent);
    }
}