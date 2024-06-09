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
#include <data/record/Water.h>

#include <random>

using namespace data::record;

template<typename Wrapped>
hydrate::RandomCulling<Wrapped>::RandomCulling::~RandomCulling() = default;

template<typename Wrapped>
std::vector<data::record::Water> hydrate::RandomCulling<Wrapped>::cull(std::vector<grid::GridMember<Water>>& placed_water) const {
    std::shuffle(placed_water.begin(), placed_water.end(), std::mt19937{std::random_device{}()}); // shuffle the molecules
    return Wrapped::cull(placed_water);
}

template class hydrate::RandomCulling<hydrate::CounterCulling>;
template class hydrate::RandomCulling<hydrate::OutlierCulling>;
template class hydrate::RandomCulling<hydrate::BodyCounterCulling>;