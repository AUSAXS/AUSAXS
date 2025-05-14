#pragma once

#include <grid/Grid.h>
#include <grid/detail/GridExcludedVolume.h>
#include <utility/observer_ptr.h>

namespace ausaxs::grid::exv {
    GridExcludedVolume create(observer_ptr<grid::Grid> grid);
}