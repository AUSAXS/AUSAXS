#pragma once

#include <vector>
#include <math/MathFwd.h>

namespace ausaxs::grid::exv {
    /**
     * @brief A simple representation of a grid-based excluded volume.
     * 
     */
    struct GridExcludedVolume {
        std::vector<Vector3<double>> interior;
        std::vector<Vector3<double>> surface;
    };
}