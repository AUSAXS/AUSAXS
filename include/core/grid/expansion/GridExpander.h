#pragma once

#include <data/atoms/AtomFF.h>
#include <data/atoms/Water.h>
#include <grid/Grid.h>
#include <grid/detail/GridMember.h>
#include <grid/expansion/SphericalExpander.h>
#include <utility/observer_ptr.h>
#include <settings/GridSettings.h>

#include <stdexcept>

namespace ausaxs::grid::volume {
    template<valid_gridmember T>
    void expand(observer_ptr<grid::Grid> grid, T& atom);

    template<valid_gridmember T>
    void deflate(observer_ptr<grid::Grid> grid, T& atom);
}

template<ausaxs::grid::valid_gridmember T>
inline void ausaxs::grid::volume::expand(observer_ptr<grid::Grid> grid, T& atom) {
    switch (settings::grid::exv::expansion_strategy) {
        case settings::grid::exv::ExvType::AtomicOnly:
            return ausaxs::grid::volume::AtomicExpander::expand_volume(grid, atom);
        case settings::grid::exv::ExvType::AtomicAndWater:
            return ausaxs::grid::volume::AtomicAndWaterExpander::expand_volume(grid, atom);
        default:
            throw std::runtime_error("GridExpander: Reached end of function. Did you forget to add a case?");
    }
}

template<ausaxs::grid::valid_gridmember T>
inline void ausaxs::grid::volume::deflate(observer_ptr<grid::Grid> grid, T& atom) {
    switch (settings::grid::exv::expansion_strategy) {
        case settings::grid::exv::ExvType::AtomicOnly:
            return ausaxs::grid::volume::AtomicExpander::deflate_volume(grid, atom);
        case settings::grid::exv::ExvType::AtomicAndWater:
            return ausaxs::grid::volume::AtomicAndWaterExpander::deflate_volume(grid, atom);
        default:
            throw std::runtime_error("GridExpander: Reached end of function. Did you forget to add a case?");
    }
}