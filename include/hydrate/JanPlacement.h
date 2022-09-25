#pragma once

#include <hydrate/PlacementStrategy.h>
#include <hydrate/Grid.h>

namespace grid {
    /**
     * @brief This strategy iterates through all bins, and for every bin which is part of the volume of an atom, it attempts to place a
     * water molecule at x±r, y±r, and z±r. If the location is valid, the molecule will be placed. This will typically generate
     * a lot of molecules, and so a culling method may be useful afterwards. 
     * 
     * The radius r is defined as the sum of @a ra and @a rh.
     */
    class JanPlacement : public PlacementStrategy {
        public:
            using PlacementStrategy::PlacementStrategy; // inherit constructor
            ~JanPlacement() override {}

            std::vector<GridMember<Hetatom>> place() const override;
        };
}