// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <gcem.hpp>

#include <cmath>
#include <numbers>
#include <version>

namespace ausaxs::gcem {
    // Fast, low-accuracy trigonometric approximations usable in constexpr contexts. Prefer these
    // only where the ~0.001 maximum error is acceptable and compile-time evaluation or raw speed
    // matters more than precision; otherwise use the constexpr_math alias defined below.
    namespace fast {
        /**
         * @brief Fast cosine function with a maximum error of 0.001.
         */
        [[maybe_unused]] inline constexpr double cos(double x) noexcept {
            // this is taken directly from https://stackoverflow.com/a/28050328
            constexpr double tp = 1./(2*std::numbers::pi);
            x *= tp;
            x -= 0.25 + ::gcem::floor(x + 0.25);
            x *= 16.*(::gcem::abs(x) - 0.5);
            x += .225*x*(::gcem::abs(x) - 1.);
            return x;
        }

        /**
         * @brief Fast sine function with a maximum error of 0.001.
         */
        inline constexpr double sin(double x) noexcept {
            // this is a slightly altered version of https://stackoverflow.com/a/28050328
            constexpr double tp = 1./(2*std::numbers::pi);
            x *= tp;
            x -= 0.5 + ::gcem::floor(x);
            x *= 16.*(::gcem::abs(x) - 0.5);
            x += .225*x*(::gcem::abs(x) - 1.);
            return x;
        }
    }
}

namespace ausaxs {
    // Alias for a constexpr-capable cmath: the standard library where it provides constexpr math
    // functions, otherwise the third-party gcem library as a fallback.
    #if defined __cpp_lib_constexpr_cmath
        namespace constexpr_math = std;
    #else
        namespace constexpr_math = ::gcem;
    #endif
}