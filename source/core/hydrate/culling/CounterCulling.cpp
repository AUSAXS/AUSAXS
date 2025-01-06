/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <hydrate/culling/CounterCulling.h>
#include <data/Molecule.h>
#include <data/atoms/Water.h>
#include <grid/detail/GridMember.h>
#include <grid/Grid.h>

#include <cmath>

using namespace ausaxs;

void hydrate::CounterCulling::cull(std::span<grid::GridMember<data::Water>>& placed_water) const {
    if (target_count == 0) {return;}
    int factor = std::floor(placed_water.size()/target_count); // reduction factor
    if (factor < 2) {return;}

    std::vector<bool> remove(placed_water.size(), false); // the water molecules which will be removed
    unsigned int pw_index = 0; // current index in placed_water
    unsigned int counter = 0; // counter
    for (int i = 0; i < static_cast<int>(placed_water.size()); ++i) {
        counter++;
        if (counter % factor != 0) {
            remove[pw_index++] = true;
            continue;
        }
    }
    molecule->get_grid()->remove_waters(remove);
}