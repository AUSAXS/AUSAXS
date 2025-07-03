// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <concepts>
#include <cstddef>

namespace ausaxs {
    template<typename C>
    concept container_type = requires(C c, int i) {
        {c[i]};
        {c.size()};
    };

    template<typename C>
    concept numeric = std::is_arithmetic_v<C>;
}