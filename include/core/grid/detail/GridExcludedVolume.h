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

        /**
         * @brief The spacing in Ångström of the cubic lattice the points lie on, or 0 if they are not lattice-supported.
         *
         * The grid-based models emit voxel centers, so both point sets are exact subsets of the sites of a single cubic
         * lattice. Recording the spacing lets consumers exploit that structure; see hist::detail::lattice, which replaces
         * the quadratic self-correlation loop with a transform on the strength of it.
         */
        double spacing = 0;

        bool has_surface() const;
        void save(const io::File& file) const;
    };
}