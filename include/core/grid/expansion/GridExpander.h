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
void ausaxs::grid::volume::expand(observer_ptr<grid::Grid> grid, T& atom) {
    switch (settings::grid::exv::expansion_strategy) {
        case settings::grid::exv::Expansion::Minimal:
            return ausaxs::grid::volume::MinimalExpander::expand_volume(grid, atom);
        case settings::grid::exv::Expansion::Full:
            return ausaxs::grid::volume::FullExpander::expand_volume(grid, atom);
        default:
            throw std::runtime_error("GridExpander: Reached end of function. Did you forget to add a case?");
    }
}

template<ausaxs::grid::valid_gridmember T>
void ausaxs::grid::volume::deflate(observer_ptr<grid::Grid> grid, T& atom) {
    switch (settings::grid::exv::expansion_strategy) {
        case settings::grid::exv::Expansion::Minimal:
            return ausaxs::grid::volume::MinimalExpander::deflate_volume(grid, atom);
        case settings::grid::exv::Expansion::Full:
            return ausaxs::grid::volume::FullExpander::deflate_volume(grid, atom);
        default:
            throw std::runtime_error("GridExpander: Reached end of function. Did you forget to add a case?");
    }
}