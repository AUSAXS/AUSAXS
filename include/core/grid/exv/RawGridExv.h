#pragma once

#include <grid/Grid.h>
#include <utility/observer_ptr.h>
#include <settings/GridSettings.h>

namespace ausaxs::grid::exv {
    struct RawGridExv {
        GridExcludedVolume create(observer_ptr<grid::Grid> grid);
    };
}

inline ausaxs::grid::exv::GridExcludedVolume ausaxs::grid::exv::RawGridExv::create(observer_ptr<grid::Grid> grid) {
    std::vector<Vector3<double>> atoms;
    atoms.reserve(grid->get_volume());

    int stride = std::round(settings::grid::exv::width/settings::grid::cell_width);
    auto[imin, imax] = grid->bounding_box_index();
    for (int i = imin.x(); i < imax.x(); i += stride) {
        for (int j = imin.y(); j < imax.y(); j += stride) {
            for (int k = imin.z(); k < imax.z(); k += stride) {
            }       
        }    
    }
}