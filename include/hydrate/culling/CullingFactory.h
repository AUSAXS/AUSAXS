#pragma once

#include <hydrate/culling/CullingStrategy.h>
#include <grid/detail/GridInternalFwd.h>
#include <settings/MoleculeSettings.h>

#include <memory>

namespace hydrate {
    namespace factory {
        /**
         * @brief Prepare a culling class. 
         */
        std::unique_ptr<CullingStrategy> construct_culling_strategy(observer_ptr<grid::Grid> grid);

        /**
         * @brief Prepare a culling class. 
         */
        std::unique_ptr<CullingStrategy> construct_culling_strategy(observer_ptr<grid::Grid> grid, const settings::hydrate::CullingStrategy& choice);
    }
}