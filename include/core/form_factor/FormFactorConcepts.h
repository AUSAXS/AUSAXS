// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <concepts>

namespace ausaxs {
    /**
     * @brief A form factor which can be evaluated at a given q value.
     */
    template<typename T>
    concept FormFactorType = requires(const T& t, double q) {
        {t.evaluate(q)} -> std::convertible_to<double>;
    };
}
