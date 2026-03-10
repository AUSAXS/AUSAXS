// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/detail/data/IntrinsicMacros.h>

namespace ausaxs::hist::detail {
    static inline float squared_dot_product(const float* v1, const float* v2) noexcept {
        float dx = v1[0] - v2[0];
        float dy = v1[1] - v2[1];
        float dz = v1[2] - v2[2];
        return dx*dx + dy*dy + dz*dz;
    }
}
