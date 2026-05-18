// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <concepts>
#include <cstddef>

namespace ausaxs {
    /**
     * @brief Requires a type to be indexable with operator[] and to expose a size().
     */
    template<typename C>
    concept container_type = requires(C c, int i) {
        {c[i]};
        {c.size()};
    };

    /**
     * @brief Requires a type to be an arithmetic (integral or floating-point) type.
     */
    template<typename C>
    concept numeric = std::is_arithmetic_v<C>;
}