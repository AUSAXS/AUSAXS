#pragma once

#include <grid/placement/PlacementStrategy.h>
#include <grid/detail/GridInternalFwd.h>
#include <settings/GridSettings.h>

#include <memory>

namespace grid {
    namespace factory {
        /**
         * @brief Prepare a placement class.
         */
        std::unique_ptr<PlacementStrategy> construct_placement_strategy(Grid* grid);

        /**
         * @brief Prepare a placement class.
         */
        std::unique_ptr<PlacementStrategy> construct_placement_strategy(Grid* grid, const settings::grid::PlacementStrategy& choice);
    }
}