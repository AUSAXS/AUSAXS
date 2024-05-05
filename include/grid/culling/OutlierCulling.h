#pragma once

#include <grid/culling/CullingStrategy.h>

namespace grid {
    /**
     * @brief Iterate through all water molecules, and count how many other molecules are nearby. Atoms counts as +1, while other water molecules counts as -2. 
     *        Then start removing the most negative water molecules until the desired count is reached. 
     */
    class OutlierCulling : public CullingStrategy {
        public:
            using CullingStrategy::CullingStrategy;
            ~OutlierCulling() override = default;

            // runs in O(n ln n) where n is the number of water molecules
            std::vector<data::record::Water> cull(std::vector<GridMember<data::record::Water>>& placed_water) const override;
        };
}