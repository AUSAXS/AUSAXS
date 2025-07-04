// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <math/MathConcepts.h>

namespace ausaxs {
    template<typename C>
    concept indexable = requires(C c, int i) {
        c[i];
    };
}