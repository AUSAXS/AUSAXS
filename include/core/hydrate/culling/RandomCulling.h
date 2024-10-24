#pragma once

#include <hydrate/culling/CullingStrategy.h>

namespace ausaxs::hydrate {
    /**
     * @brief Iterates through all placed water molecules, rejecting all but the nth, where n is determined from the desired number of water molecules. 
     */
    template<typename Wrapped>
    class RandomCulling : public Wrapped {
        public:
            using Wrapped::Wrapped;
            ~RandomCulling() override;

            // runs in O(n) where n is the number of water molecules
            std::vector<data::record::Water> cull(std::vector<grid::GridMember<data::record::Water>>& placed_water) const override;
    };
}