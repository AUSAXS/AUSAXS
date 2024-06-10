#pragma once

#include <math/MathConcepts.h>

template<typename C>
concept indexable = requires(C c, int i) {
    c[i];
};