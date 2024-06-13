/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <hydrate/culling/NoCulling.h>
#include <grid/detail/GridMember.h>
#include <data/record/Water.h>
#include <grid/Grid.h>

using namespace data::record;

std::vector<data::record::Water> hydrate::NoCulling::cull(std::vector<grid::GridMember<Water>>& placed_water) const {
    std::vector<Water> final_water(placed_water.size());
    std::transform(placed_water.begin(), placed_water.end(), final_water.begin(), [] (const grid::GridMember<Water>& gm) -> const auto& {return gm.get_atom();});
    return final_water;
}