#include <grid/exv/GridExvStrategy.h>
#include <grid/exv/RawGridExv.h>
#include <grid/exv/RawGridWithSurfaceExv.h>
#include <utility/Logging.h>
#include <settings/ExvSettings.h>

using namespace ausaxs;
using namespace ausaxs::grid::exv;

GridExcludedVolume grid::exv::create(observer_ptr<grid::Grid> grid) {
    switch (settings::exv::exv_method) {
        case settings::exv::ExvMethod::GridSurface:
            return RawGridWithSurfaceExv::create(grid);
        case settings::exv::ExvMethod::Grid:
        case settings::exv::ExvMethod::GridScalable:
        case settings::exv::ExvMethod::WAXSiS:
            return RawGridExv::create(grid);

        default:
            logging::log("GridExvStrategy::create: Chosen exv model does not use a grid-based excluded volume. Returning empty object.");
            return {{}, {}};
    }
}