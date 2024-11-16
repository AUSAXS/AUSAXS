#pragma once

#include <hydrate/culling/CullingStrategy.h>

namespace ausaxs::hydrate {
    /**
     * @brief No culling will be performed.  
     */
    class NoCulling : public CullingStrategy {
        public:
            using CullingStrategy::CullingStrategy;
            ~NoCulling() override = default;

            std::vector<data::record::Water> cull(std::vector<grid::GridMember<data::record::Water>>& placed_water) const override;
    };       
}