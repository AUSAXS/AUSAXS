// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <gcem.hpp>

#include <version>

namespace ausaxs {
    // Alias for a constexpr-capable cmath: the standard library where it provides constexpr math
    // functions, otherwise the third-party gcem library as a fallback.
    #if defined __cpp_lib_constexpr_cmath
        namespace constexpr_math = std;
    #else
        namespace constexpr_math = ::gcem;
    #endif
}