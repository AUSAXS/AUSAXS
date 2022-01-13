#include "hydrate/CullingStrategy.h"
#include "hydrate/Grid.h"

/**
 * @brief Iterates through all placed water molecules, rejecting all but the nth, where n is determined from the desired number of water molecules. 
 */
class CounterCulling : public CullingStrategy {
public:
    using CullingStrategy::CullingStrategy;
    ~CounterCulling() override {}

    // runs in O(n) where n is the number of water molecules
    vector<Hetatom> cull(vector<GridMember<Hetatom>>& placed_water) const override {
        std::cout << "crash? cull1" << std::endl;
        if (target_count == 0) {
            vector<Hetatom> final_water(placed_water.size());
            std::transform(placed_water.begin(), placed_water.end(), final_water.begin(), [] (GridMember<Hetatom>& gm) {return gm.atom;});
            return final_water;
        }

        std::cout << "crash? cull2" << std::endl;
        int factor = std::floor(placed_water.size()/target_count); // reduction factor
        if (factor < 2) {
            vector<Hetatom> final_water(placed_water.size());
            std::transform(placed_water.begin(), placed_water.end(), final_water.begin(), [] (GridMember<Hetatom>& gm) {return gm.atom;});
            return final_water;
        }

        std::cout << "crash? cull3" << std::endl;
        vector<Hetatom> final_water(placed_water.size()); // the final water molecules that will be used
        vector<Hetatom> removed_water(placed_water.size()); // the water molecules which will be removed
        size_t rm_index = 0; // current index in removed_water
        size_t pw_index = 0; // current index in placed_water
        size_t counter = 0; // counter
        std::cout << "crash? cull4" << std::endl;
        for (const auto& a : placed_water) {
            counter++;
            if (counter % factor != 0) {
                removed_water[rm_index++] = a.atom;
                continue;
            }
            final_water[pw_index++] = a.atom;
        }
        std::cout << "crash? cull5" << std::endl;
        removed_water.resize(rm_index);
        final_water.resize(pw_index);
        std::cout << "crash? cull6" << std::endl;
        grid->remove(removed_water);
        std::cout << "crash? cull7" << std::endl;
        return final_water;
    }
};