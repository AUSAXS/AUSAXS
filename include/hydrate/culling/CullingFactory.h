#pragma once

#include <hydrate/culling/CullingStrategy.h>
#include <utility/SmartOption.h>

namespace settings {
    namespace grid {
        enum class CullingStrategy {CounterStrategy, OutlierStrategy, RandomStrategy};

        extern settings::detail::SmartOption<CullingStrategy> culling_strategy;
    }
}

namespace grid {
    namespace factory {
        /**
         * @brief Prepare a culling class. 
         */
        std::unique_ptr<CullingStrategy> construct_culling_strategy(Grid* grid, settings::grid::CullingStrategy choice = settings::grid::culling_strategy);
    }
}