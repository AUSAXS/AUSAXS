/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <hydrate/culling/CounterCulling.h>
#include <grid/detail/GridMember.h>
#include <data/record/Water.h>
#include <grid/Grid.h>

#include <cmath>

using namespace data::record;

std::vector<data::record::Water> hydrate::CounterCulling::cull(std::vector<Water>&& placed_water) const {
    if (target_count == 0) {
        return placed_water;
    }

    int factor = std::floor(placed_water.size()/target_count); // reduction factor
    if (factor < 2) {
        return placed_water;
    }

    std::vector<Water> final_water; // the final water molecules that will be used
    std::vector<Water> removed_water; // the water molecules which will be removed
    final_water.reserve(placed_water.size());
    removed_water.reserve(placed_water.size());
    unsigned int counter = 0;
    for (const auto& a : placed_water) {
        counter++;
        if (counter % factor != 0) {
            removed_water.emplace_back(a);
            continue;
        }
        final_water.emplace_back(a);
    }
    grid->add(final_water);
    return final_water;
}