#pragma once

#include <grid/expansion/GridExpander.h>
#include <settings/GridSettings.h>

#include <memory>

namespace ausaxs::grid::expander::factory {
    std::unique_ptr<GridExpander> construct(observer_ptr<Grid> grid);
    std::unique_ptr<GridExpander> construct(observer_ptr<Grid> grid, settings::grid::exv::Expansion choice);
}