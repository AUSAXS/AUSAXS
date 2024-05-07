/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <hydrate/culling/RandomCulling.h>
#include <grid/detail/GridMember.h>
#include <grid/Grid.h>
#include <data/record/Water.h>

#include <random>

using namespace data::record;

hydrate::RandomCulling::~RandomCulling() = default;

std::vector<data::record::Water> hydrate::RandomCulling::cull(std::vector<Water>&& placed_water) const {
    std::shuffle(placed_water.begin(), placed_water.end(), std::mt19937{std::random_device{}()}); // shuffle the molecules
    grid->add(placed_water);
    return CounterCulling::cull(std::move(placed_water));
}
