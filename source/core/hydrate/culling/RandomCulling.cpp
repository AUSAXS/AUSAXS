// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

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
void hydrate::RandomCulling<Wrapped>::cull(std::span<grid::GridMember<data::Water>>& placed_water) const {
    std::shuffle(placed_water.begin(), placed_water.end(), std::mt19937{std::random_device{}()}); // shuffle the molecules
    Wrapped::cull(placed_water);
}

template class ausaxs::hydrate::RandomCulling<hydrate::CounterCulling>;
template class ausaxs::hydrate::RandomCulling<hydrate::OutlierCulling>;
template class ausaxs::hydrate::RandomCulling<hydrate::BodyCounterCulling>;