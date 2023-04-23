#pragma once

#include <hydrate/placement/PlacementStrategy.h>
#include <hydrate/GridSettings.h>

namespace grid {
    namespace factory {
        /**
         * @brief Prepare a placement class.
         */
        std::unique_ptr<PlacementStrategy> construct_placement_strategy(Grid* grid, settings::grid::PlacementStrategy choice = settings::grid::placement_strategy);
    }
}