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