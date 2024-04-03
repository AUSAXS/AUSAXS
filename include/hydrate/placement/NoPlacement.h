#pragma once

#include <hydrate/placement/PlacementStrategy.h>

namespace grid {
    /**
     * @brief This strategy will not place any water molecules. 
     */
    class NoPlacement : public PlacementStrategy {
        public:
            using PlacementStrategy::PlacementStrategy; // inherit constructor
            ~NoPlacement() override = default;

            std::vector<GridMember<data::record::Water>> place() const override;
        };
}