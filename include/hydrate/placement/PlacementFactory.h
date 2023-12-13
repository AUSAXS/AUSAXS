#pragma once

#include <hydrate/placement/PlacementStrategy.h>
#include <hydrate/detail/GridInternalFwd.h>

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