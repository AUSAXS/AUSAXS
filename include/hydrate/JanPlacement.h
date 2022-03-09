#include "hydrate/PlacementStrategy.h"
#include "hydrate/Grid.h"

namespace grid {
    /**
     * @brief JanStrategy
     * 
     * This strategy iterates through all bins, and for every bin which is part of the volume of an atom, it attempts to place a
     * water molecule at x±r, y±r, and z±r. If the location is valid, the molecule will be placed. This will typically generate
     * a lot of molecules, and so a culling method may be useful afterwards. 
     * 
     * It only uses the atomic radius @a ra, and ignores the value set to the hydration atoms @a ra. 
     * 
     */
    class JanPlacement : public PlacementStrategy {
    public:
        using PlacementStrategy::PlacementStrategy; // inherit constructor
        ~JanPlacement() override {}

        vector<GridMember<Hetatom>> place() const override;
    };
}