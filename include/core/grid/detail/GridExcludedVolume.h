// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <math/MathFwd.h>
#include <io/IOFwd.h>

#include <vector>

namespace ausaxs::grid::exv {
    /**
     * @brief A simple representation of a grid-based excluded volume.
     * 
     */
    struct GridExcludedVolume {
        std::vector<Vector3<double>> interior;
        std::vector<Vector3<double>> surface;

        bool has_surface() const;
        void save(const io::File& file) const;
    };
}