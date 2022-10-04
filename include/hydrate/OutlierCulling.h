#pragma once

#include <hydrate/CullingStrategy.h>
#include <hydrate/Grid.h>

namespace grid {
    /**
     * @brief Iterate through all water molecules, and count how many other molecules are nearby. Atoms counts as +1, while other water molecules counts as -2. 
     *        Then start removing the most negative water molecules until the desired count is reached. 
     */
    class OutlierCulling : public CullingStrategy {
        public:
            using CullingStrategy::CullingStrategy;
            ~OutlierCulling() override {}

            // runs in O(n ln n) where n is the number of water molecules
            std::vector<Water> cull(std::vector<GridMember<Water>>& placed_water) const override;
        };
}