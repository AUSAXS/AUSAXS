// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hydrate/generation/HydrationFactory.h>
#include <hydrate/generation/JanHydration.h>
#include <hydrate/generation/RadialHydration.h>
#include <hydrate/generation/AxesHydration.h>
#include <hydrate/generation/PepsiHydration.h>
#include <hydrate/generation/NoHydration.h>
#include <hydrate/culling/CullingFactory.h>
#include <settings/MoleculeSettings.h>
#include <utility/Exceptions.h>
#include <math/Vector3.h>

using namespace ausaxs;
using namespace ausaxs::hydrate;

std::unique_ptr<HydrationStrategy> factory::construct_hydration_generator(observer_ptr<data::Molecule> protein) {
    return factory::construct_hydration_generator(protein, settings::hydrate::hydration_strategy, settings::hydrate::culling_strategy);
}

std::unique_ptr<HydrationStrategy> factory::construct_hydration_generator(observer_ptr<data::Molecule> protein, settings::hydrate::HydrationStrategy choice) {
    return factory::construct_hydration_generator(protein, choice, settings::hydrate::culling_strategy);
}

std::unique_ptr<HydrationStrategy> factory::construct_hydration_generator(
    observer_ptr<data::Molecule> protein, settings::hydrate::HydrationStrategy choice, settings::hydrate::CullingStrategy culling_strategy) 
{
    switch (choice) {
        case settings::hydrate::HydrationStrategy::AxesStrategy: 
            return std::make_unique<AxesHydration>(protein, factory::construct_culling_strategy(protein, culling_strategy));
        case settings::hydrate::HydrationStrategy::RadialStrategy:
            return std::make_unique<RadialHydration>(protein, factory::construct_culling_strategy(protein, culling_strategy));
        case settings::hydrate::HydrationStrategy::JanStrategy: 
            return std::make_unique<JanHydration>(protein, factory::construct_culling_strategy(protein, culling_strategy));
        case settings::hydrate::HydrationStrategy::NoStrategy: 
            return std::make_unique<NoHydration>();
        case settings::hydrate::HydrationStrategy::PepsiStrategy:
            return std::make_unique<PepsiHydration>(protein, factory::construct_culling_strategy(protein, culling_strategy));
        default: 
            throw except::unknown_argument("hydrate::factory::construct_hydration_generator: Unkown HydrationStrategy. Did you forget to add it to the switch statement?");
    }
}