#pragma once

#include <math/MathConcepts.h>

namespace ausaxs {
    template<typename C>
    concept indexable = requires(C c, int i) {
        c[i];
    };
}