#pragma once

#include <hydrate/PlacementStrategy.h>
#include <hydrate/Grid.h>

namespace grid {
    /**
     * @brief This strategy iterates through all atoms and attempts to place a water molecule at x±r, y±r, and z±r. 
     *        If the location is valid, the molecule will be placed. This will typically generate a lot of molecules, 
     *        and so a culling method may be useful afterwards. 
     */
    class AxesPlacement : public PlacementStrategy {
    public:
        using PlacementStrategy::PlacementStrategy; // inherit constructor
        ~AxesPlacement() override {}

        vector<GridMember<Hetatom>> place() const override;

    private:

        /**
         * @brief Check if a water molecule can be placed at the given location. 
         * @param loc the location to be checked. 
         * @return True if this is an acceptable location, false otherwise.
         */
        bool collision_check(const vector<int> loc) const;
    };
}