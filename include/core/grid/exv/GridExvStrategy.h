// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <grid/Grid.h>
#include <grid/detail/GridExcludedVolume.h>
#include <utility/observer_ptr.h>

namespace ausaxs::grid::exv {
    GridExcludedVolume create(observer_ptr<grid::Grid> grid);

    /**
     * @brief The stride in grid cells between neighbouring excluded volume points.
     */
    int point_stride();

    /**
     * @brief The spacing in Ångström of the cubic lattice all grid-based excluded volume points lie on.
     */
    double point_spacing();
}