#pragma once

#include <grid/Grid.h>
#include <grid/detail/GridExcludedVolume.h>
#include <utility/observer_ptr.h>

namespace ausaxs::grid::exv {
    GridExcludedVolume create(observer_ptr<grid::Grid> grid);
}

inline ausaxs::grid::exv::GridExcludedVolume ausaxs::grid::exv::create(observer_ptr<grid::Grid> grid) {
    
}