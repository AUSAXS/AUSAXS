#pragma once

#include <hydrate/generation/GridBasedHydration.h>
#include <math/MathFwd.h>

#include <vector>

namespace hydrate {
    /**
     * @brief This strategy iterates through all atoms, and attempts to place a water molecule at a distance r along a number of radial lines originating from each atom. 
     * For each possible location, a suitability score is calculated, which favors the surface of the molecule, and penalizes cavities (including internal spaces). 
     * 
     * Although more calculations are involved for each location than the AxesPlacement strategy, the complexity is the same.
     * 
     * The radius r is defined as the sum of @a ra and @a rh.
     */
    class RadialHydration : public GridBasedHydration {
        public:
            RadialHydration(observer_ptr<data::Molecule> protein);
            virtual ~RadialHydration();

            std::vector<data::record::Water> generate_explicit_hydration() override;

        private:
            observer_ptr<grid::Grid> grid;

            void initialize() override;
            void prepare_rotations(int divisions = 8);

            std::vector<Vector3<int>> rot_bins_1rh; // rotation bins at 1rh radius
            std::vector<Vector3<int>> rot_bins_3rh; // rotation bins at 3rh radius
            std::vector<Vector3<int>> rot_bins_5rh; // rotation bins at 5rh radius
            std::vector<Vector3<int>> rot_bins_7rh; // rotation bins at 7rh radius
            std::vector<Vector3<double>> rot_locs;  // absolute locations of the rotation bins

            /**
             * @brief Check if a water molecule can be placed at the given location. 
             * @param loc the location to be checked. 
             * @param skip_bin location to be excluded from the check. 
             * @return True if this is an acceptable location, false otherwise.
             */
            bool collision_check(const Vector3<int>& loc, const Vector3<int>& skip_bin) const;
    };
}