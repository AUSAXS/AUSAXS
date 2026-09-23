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
        std::vector<Vector3<double>> interior;
        std::vector<Vector3<double>> surface;

        /**
         * @brief The integer lattice sites of @a interior and @a surface, index for index.
         *
         * The grid-based models emit voxel centers, so both point sets are exact subsets of the sites of a single cubic
         * lattice. The sites let consumers exploit that structure directly; see hist::detail::lattice, which replaces
         * the quadratic self-correlation loop with a transform on the strength of it.
         */
        std::vector<Vector3<int>> interior_sites;
        std::vector<Vector3<int>> surface_sites;

        // the spacing in Ångström of the lattice the sites are expressed on
        double spacing = 0;

        bool has_surface() const;
        void save(const io::File& file) const;
    };
}