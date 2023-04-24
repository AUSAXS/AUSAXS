#pragma once

#include <hydrate/culling/CullingStrategy.h>
#include <settings/GridSettings.h>

namespace grid {
    namespace factory {
        /**
         * @brief Prepare a culling class. 
         */
        std::unique_ptr<CullingStrategy> construct_culling_strategy(Grid* grid, settings::grid::CullingStrategy choice = settings::grid::culling_strategy);
    }
}