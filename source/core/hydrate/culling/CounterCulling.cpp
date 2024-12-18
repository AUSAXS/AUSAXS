/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <hydrate/culling/CounterCulling.h>
#include <data/Molecule.h>
#include <grid/detail/GridMember.h>
#include <grid/Grid.h>

#include <cmath>

using namespace ausaxs;

std::vector<grid::GridMember<data::Water>>& hydrate::CounterCulling::cull(std::vector<grid::GridMember<data::Water>>& placed_water) const {
    auto return_input = [&placed_water] () {
        std::vector<data::Water> final_water(placed_water.size());
        std::transform(placed_water.begin(), placed_water.end(), final_water.begin(), [] (const grid::GridMember<data::Water>& gm) -> const auto& {return gm.get_atom();});
        return final_water;
    };

    if (target_count == 0) {return molecule->get_grid()->w_members;}
    int factor = std::floor(placed_water.size()/target_count); // reduction factor
    if (factor < 2) {return molecule->get_grid()->w_members;}

    std::vector<bool> remove(placed_water.size(), false); // the water molecules which will be removed
    unsigned int rm_index = 0; // current index in removed_water
    unsigned int pw_index = 0; // current index in placed_water
    unsigned int counter = 0; // counter
    for (int i = 0; i < static_cast<int>(placed_water.size()); ++i) {
        counter++;
        if (counter % factor != 0) {
            remove[pw_index++] = true;
            continue;
        }
    }
    molecule->get_grid()->remove(removed_water);
    return final_water;
}