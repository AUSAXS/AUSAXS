#include <grid/exv/GridExvStrategy.h>
#include <grid/exv/RawGridExv.h>
#include <grid/exv/RawGridWithSurfaceExv.h>
#include <settings/HistogramSettings.h>

using namespace ausaxs;
using namespace ausaxs::grid::exv;

GridExcludedVolume grid::exv::create(observer_ptr<grid::Grid> grid) {
    switch (settings::hist::get_exv_strategy()) {
        case settings::hist::ExvMethod::Grid:
            return RawGridExv::create(grid);
        case settings::hist::ExvMethod::GridWithSurface:
            return RawGridWithSurfaceExv::create(grid);
        default:
            throw except::unexpected("GridExvStrategy::create: Unknown ExvMethod. Did you forget to add it to the switch statement?");
    }
}