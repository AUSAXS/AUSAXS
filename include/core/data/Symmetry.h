#pragma once

#include <math/Vector3.h>

namespace ausaxs::data::detail {
    struct Symmetry {
        // translational vector with respect to the original body
        Vector3<double> translate;

        // rotational vector with respect to the original body. Each component describes the rotation around the corresponding axis
        Vector3<double> rotate;

        // the number of times the symmetry should be repeated
        int repeat = 1;
    };
}