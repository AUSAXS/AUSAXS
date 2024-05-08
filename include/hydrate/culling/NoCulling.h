#pragma once

#include <hydrate/culling/CullingStrategy.h>

namespace hydrate {
    /**
     * @brief No culling will be performed.  
     */
    class NoCulling : public CullingStrategy {
        public:
            using CullingStrategy::CullingStrategy;
            ~NoCulling() override = default;

            std::vector<data::record::Water> cull(std::vector<data::record::Water>&& placed_water) const override;
    };       
}