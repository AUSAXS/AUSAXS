#pragma once

#include <grid/detail/GridExcludedVolume.h>
#include <grid/GridFwd.h>
#include <utility/observer_ptr.h>

namespace ausaxs::grid::exv {
    struct RawGridWithSurfaceExv {
        static GridExcludedVolume create(observer_ptr<grid::Grid> grid, bool detect_surface = true);
    };
}