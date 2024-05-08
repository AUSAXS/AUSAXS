/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <hydrate/culling/NoCulling.h>
#include <grid/detail/GridMember.h>
#include <data/record/Water.h>
#include <grid/Grid.h>

using namespace data::record;

std::vector<data::record::Water> hydrate::NoCulling::cull(std::vector<Water>&& placed_water) const {
    return placed_water;
}