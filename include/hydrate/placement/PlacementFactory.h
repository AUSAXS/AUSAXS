#pragma once

#include <hydrate/placement/PlacementStrategy.h>

#include <memory>

namespace settings::grid {enum class PlacementStrategy;}
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