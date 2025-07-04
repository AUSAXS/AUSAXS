// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hydrate/generation/GridBasedHydration.h>
#include <math/MathFwd.h>

namespace ausaxs::hydrate {
    /**
     * @brief This strategy iterates through all atoms and attempts to place a water molecule at x±r, y±r, and z±r. 
     * If the location is valid, the molecule will be placed. This will typically generate a lot of molecules, 
     * and so a culling method may be useful afterwards. 
     * 
     * The radius r is defined as the sum of @a ra and @a rh.
     */
    class AxesHydration : public GridBasedHydration {
        public:
            AxesHydration(observer_ptr<data::Molecule> protein);
            AxesHydration(observer_ptr<data::Molecule> protein, std::unique_ptr<CullingStrategy> culling_strategy);
            virtual ~AxesHydration();

            std::span<grid::GridMember<data::Water>> generate_explicit_hydration(std::span<grid::GridMember<data::AtomFF>> atoms) override;

            bool global() const override {return false;}

        private:
            void initialize() override;

            /**
            * @brief Check if a water molecule can be placed at the given location. 
            * @param loc the location to be checked. 
            * @return True if this is an acceptable location, false otherwise.
            */
            bool collision_check(const Vector3<unsigned int>& loc, double ra) const;
    };
}