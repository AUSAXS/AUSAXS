#pragma once

#include "hydrate/PlacementStrategy.h"
#include "hydrate/Grid.h"
#include "settings.h"

namespace grid {
    /**
     * @brief Description
     */
    class RadialPlacement : public PlacementStrategy {
    public:
        RadialPlacement(Grid* grid) : PlacementStrategy(grid) {
            prepare_rotations();
        }
        ~RadialPlacement() override {}

        void prepare_rotations(const int divisions = 8);

        // the vectors representing the bin offsets of rotations
        vector<vector<int>> rot_bins_1rh; // rotation bins at 1rh radius
        vector<vector<int>> rot_bins_3rh; // rotation bins at 3rh radius
        vector<vector<int>> rot_bins_5rh; // rotation bins at 5rh radius
        vector<vector<int>> rot_bins_7rh; // rotation bins at 7rh radius
        vector<vector<int>> rot_bins_rarh; // rotation bins at rarh radius
        vector<Vector3> rot_locs_rarh; // exact locations of the rarh bins

        vector<GridMember<Hetatom>> place() const override;

    private:
        /**
         * @brief Check if a water molecule can be placed at the given location. 
         * @param loc the location to be checked. 
         * @param skip_bin location to be excluded from the check. 
         * @return True if this is an acceptable location, false otherwise.
         */
        bool collision_check(const vector<int> loc, const vector<int> skip_bin) const;
    };
}