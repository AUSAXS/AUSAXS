// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <grid/Grid.h>

namespace ausaxs::em::grid {
    /**
     * @brief A customized grid for EM maps. 
     *
     * This grid has been modified to use the EM map resolution as the minimum atomic radius.
     * This is necessary to circumvent the expanded radius becoming 0 for the map atoms with the "unknown" form factor.
     */
    class EMGrid : public ausaxs::grid::Grid {
        public:
            using Grid::Grid;
            ~EMGrid() = default;

            double get_atomic_radius(form_factor::form_factor_t) const override;
    };
}