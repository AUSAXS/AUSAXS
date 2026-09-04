// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <grid/exv/GridExvStrategy.h>
#include <grid/exv/RawGridExv.h>
#include <grid/exv/RawGridWithSurfaceExv.h>
#include <utility/Logging.h>
#include <settings/ExvSettings.h>
#include <settings/GridSettings.h>

#include <algorithm>
#include <cmath>

using namespace ausaxs;
using namespace ausaxs::grid::exv;

int grid::exv::point_stride() {
    return std::max(1., std::round(settings::grid::exv::width/settings::grid::cell_width));
}

double grid::exv::point_spacing() {
    return point_stride()*settings::grid::cell_width;
}

GridExcludedVolume grid::exv::create(observer_ptr<grid::Grid> grid) {
    switch (settings::exv::exv_method) {
        case settings::exv::ExvMethod::GridSurface:
            return RawGridWithSurfaceExv::create(grid);
        case settings::exv::ExvMethod::Grid:
        case settings::exv::ExvMethod::GridScalable:
        case settings::exv::ExvMethod::WAXSiS:
            return RawGridExv::create(grid);

        default:
            logging::log("GridExvStrategy::create: Chosen exv model does not use a grid-based excluded volume. Returning empty object.");
            return {{}, {}};
    }
}