#pragma once

#include <hydrate/culling/CounterCulling.h>

namespace grid {
    /**
     * @brief Iterates through all placed water molecules, rejecting all but the nth, where n is determined from the desired number of water molecules. 
     */
    class RandomCulling : public CounterCulling {
        public:
            using CounterCulling::CounterCulling;
            ~RandomCulling() override;

            // runs in O(n) where n is the number of water molecules
            std::vector<data::record::Water> cull(std::vector<GridMember<data::record::Water>>& placed_water) const override;
    };       
}