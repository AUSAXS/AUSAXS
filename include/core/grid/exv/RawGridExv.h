// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <grid/Grid.h>
#include <grid/detail/GridObj.h>
#include <utility/observer_ptr.h>
#include <settings/GridSettings.h>

namespace ausaxs::grid::exv {
    struct RawGridExv {
        static GridExcludedVolume create(observer_ptr<grid::Grid> grid);
    };
}