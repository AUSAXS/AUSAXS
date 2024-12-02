#pragma once

#include <math/Vector3.h>

namespace ausaxs::data::detail {
    struct Symmetry {
        Vector3<double> translate;
        Vector3<double> rotate;
    };
}