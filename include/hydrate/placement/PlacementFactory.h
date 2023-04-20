#pragma once

#include <hydrate/placement/PlacementStrategy.h>
#include <utility/SmartOption.h>

namespace settings {
    namespace grid {
        enum class PlacementStrategy {
            AxesStrategy, 
            RadialStrategy, 
            JanStrategy
        };
        extern settings::detail::SmartOption<PlacementStrategy> placement_strategy;
    }
}

namespace grid {
    namespace factory {
        /**
         * @brief Prepare a placement class.
         */
        std::unique_ptr<PlacementStrategy> construct_placement_strategy(Grid* grid, settings::grid::PlacementStrategy choice = settings::grid::placement_strategy);
    }
}