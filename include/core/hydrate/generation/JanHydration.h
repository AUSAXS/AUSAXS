#pragma once

#include <hydrate/generation/GridBasedHydration.h>
#include <math/MathFwd.h>

namespace hydrate {
    /**
     * @brief This strategy iterates through all bins, and for every bin which is part of the volume of an atom, it attempts to place a
     * water molecule at x±r, y±r, and z±r. If the location is valid, the molecule will be placed. This will typically generate
     * a lot of molecules, and so a culling method may be useful afterwards. 
     * 
     * The radius r is defined as the sum of @a ra and @a rh.
     */
    class JanHydration : public GridBasedHydration {
        public:
            JanHydration(observer_ptr<data::Molecule> protein);
            JanHydration(observer_ptr<data::Molecule> protein, std::unique_ptr<CullingStrategy> culling_strategy);
            ~JanHydration() override = default;

            std::vector<grid::GridMember<data::record::Water>> generate_explicit_hydration() override;
        };
}