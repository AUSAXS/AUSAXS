// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <grid/Grid.h>
#include <grid/detail/GridObj.h>
#include <settings/GridSettings.h>
#include <utility/observer_ptr.h>

namespace ausaxs::grid::exv {
    struct RawGridExv {
        static GridExcludedVolume create(observer_ptr<grid::Grid> grid);
    };
}