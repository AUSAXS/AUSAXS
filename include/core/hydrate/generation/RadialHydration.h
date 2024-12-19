#pragma once

#include <hydrate/generation/GridBasedHydration.h>
#include <grid/detail/GridInternalFwd.h>
#include <math/MathFwd.h>

#include <vector>
#include <functional>

namespace ausaxs::hydrate {
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
            RadialHydration(observer_ptr<data::Molecule> protein, std::unique_ptr<CullingStrategy> culling_strategy);
            virtual ~RadialHydration();

            std::vector<grid::GridMember<data::Water>> generate_explicit_hydration() override;

            static void set_noise_generator(std::function<Vector3<double>()>&& noise_function);

        private:
            std::vector<Vector3<int>> rot_bins_1rh; // rotation bins at 1rh radius
            std::vector<Vector3<int>> rot_bins_3rh; // rotation bins at 3rh radius
            std::vector<Vector3<int>> rot_bins_5rh; // rotation bins at 5rh radius
            std::vector<Vector3<int>> rot_bins_7rh; // rotation bins at 7rh radius
            std::vector<Vector3<double>> rot_locs;  // absolute locations of the rotation bins
            static std::function<Vector3<double>()> noise_generator;

            void initialize() override;
            void prepare_rotations(int divisions = 8);

            /**
             * @brief Check if a water molecule can be placed at the given location. 
             *
             * @param loc the location to be checked. 
             *
             * @return True if this is an acceptable location, false otherwise.
             */
            bool collision_check(const Vector3<int>& loc) const;
    };
}