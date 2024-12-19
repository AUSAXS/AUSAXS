/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <hydrate/culling/RandomCulling.h>
#include <hydrate/culling/CounterCulling.h>
#include <hydrate/culling/OutlierCulling.h>
#include <hydrate/culling/BodyCounterCulling.h>
#include <grid/detail/GridMember.h>
#include <grid/Grid.h>

#include <random>

using namespace ausaxs;

template<typename Wrapped>
hydrate::RandomCulling<Wrapped>::RandomCulling::~RandomCulling() = default;

template<typename Wrapped>
void hydrate::RandomCulling<Wrapped>::cull(std::vector<grid::GridMember<data::Water>>& placed_water) const {
    std::shuffle(placed_water.begin(), placed_water.end(), std::mt19937{std::random_device{}()}); // shuffle the molecules
    Wrapped::cull(placed_water);
}

template class ausaxs::hydrate::RandomCulling<hydrate::CounterCulling>;
template class ausaxs::hydrate::RandomCulling<hydrate::OutlierCulling>;
template class ausaxs::hydrate::RandomCulling<hydrate::BodyCounterCulling>;