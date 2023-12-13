#pragma once

#include <hydrate/culling/CullingStrategy.h>
#include <hydrate/detail/GridInternalFwd.h>
#include <settings/GridSettings.h>

#include <memory>

namespace grid {
    namespace factory {
        /**
         * @brief Prepare a culling class. 
         */
        std::unique_ptr<CullingStrategy> construct_culling_strategy(Grid* grid);

        /**
         * @brief Prepare a culling class. 
         */
        std::unique_ptr<CullingStrategy> construct_culling_strategy(Grid* grid, const settings::grid::CullingStrategy& choice);
    }
}