#include <grid/exv/GridExvStrategy.h>
#include <grid/exv/RawGridExv.h>
#include <grid/exv/RawGridWithSurfaceExv.h>
#include <settings/GridSettings.h>

using namespace ausaxs;
using namespace ausaxs::grid::exv;

GridExcludedVolume grid::exv::create(observer_ptr<grid::Grid> grid) {
    switch (settings::grid::exv::exv_strategy) {
        case settings::grid::exv::Exv::Grid:
            return RawGridExv::create(grid);
        case settings::grid::exv::Exv::GridWithSurface:
            return RawGridWithSurfaceExv::create(grid);
        default:
            throw except::io_error("Invalid grid excluded volume strategy.");
    }
}