/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <grid/culling/RandomCulling.h>
#include <grid/GridMember.h>
#include <data/record/Water.h>

#include <random>

using namespace grid;
using namespace data::record;

RandomCulling::~RandomCulling() = default;

std::vector<data::record::Water> RandomCulling::cull(std::vector<GridMember<Water>>& placed_water) const {
    std::shuffle(placed_water.begin(), placed_water.end(), std::mt19937{std::random_device{}()}); // shuffle the molecules
    return CounterCulling::cull(placed_water);
}
