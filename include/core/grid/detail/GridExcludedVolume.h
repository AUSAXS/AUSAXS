// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <io/IOFwd.h>
#include <math/MathFwd.h>

#include <vector>

namespace ausaxs::grid::exv {
    /**
     * @brief A simple representation of a grid-based excluded volume.
     * 
     */
    struct GridExcludedVolume {
        std::vector<Vector3<double>> interior;    // real interior positions
        std::vector<Vector3<double>> surface;     // real surface positions
        std::vector<Vector3<int>> interior_sites; // interior indices
        std::vector<Vector3<int>> surface_sites;  // surface indices
        double spacing = 0; // cell spacing of the interior/surface indices

        bool has_surface() const;
        void save(const io::File& file) const;
    };
}