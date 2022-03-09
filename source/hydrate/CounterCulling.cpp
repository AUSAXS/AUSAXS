#include "hydrate/CounterCulling.h"
#include "hydrate/Grid.h"

vector<Hetatom> grid::CounterCulling::cull(vector<grid::GridMember<Hetatom>>& placed_water) const {
    if (target_count == 0) {
        vector<Hetatom> final_water(placed_water.size());
        std::transform(placed_water.begin(), placed_water.end(), final_water.begin(), [] (GridMember<Hetatom>& gm) {return gm.atom;});
        return final_water;
    }

    int factor = std::floor(placed_water.size()/target_count); // reduction factor
    if (factor < 2) {
        vector<Hetatom> final_water(placed_water.size());
        std::transform(placed_water.begin(), placed_water.end(), final_water.begin(), [] (GridMember<Hetatom>& gm) {return gm.atom;});
        return final_water;
    }

    vector<Hetatom> final_water(placed_water.size()); // the final water molecules that will be used
    vector<Hetatom> removed_water(placed_water.size()); // the water molecules which will be removed
    size_t rm_index = 0; // current index in removed_water
    size_t pw_index = 0; // current index in placed_water
    size_t counter = 0; // counter
    for (const auto& a : placed_water) {
        counter++;
        if (counter % factor != 0) {
            removed_water[rm_index++] = a.atom;
            continue;
        }
        final_water[pw_index++] = a.atom;
    }
    removed_water.resize(rm_index);
    final_water.resize(pw_index);
    grid->remove(removed_water);
    return final_water;
}