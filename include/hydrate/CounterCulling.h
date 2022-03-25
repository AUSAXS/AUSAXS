#pragma once

#include "hydrate/CullingStrategy.h"
#include "hydrate/Grid.h"

namespace grid {
    /**
     * @brief Iterates through all placed water molecules, rejecting all but the nth, where n is determined from the desired number of water molecules. 
     */
    class CounterCulling : public CullingStrategy {
    public:
        using CullingStrategy::CullingStrategy;
        ~CounterCulling() override {}

        // runs in O(n) where n is the number of water molecules
        vector<Hetatom> cull(vector<GridMember<Hetatom>>& placed_water) const override;
    };       
}