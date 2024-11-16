/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <hydrate/culling/CullingFactory.h>
#include <hydrate/culling/CounterCulling.h>
#include <hydrate/culling/BodyCounterCulling.h>
#include <hydrate/culling/OutlierCulling.h>
#include <hydrate/culling/RandomCulling.h>
#include <hydrate/culling/NoCulling.h>
#include <settings/GridSettings.h>
#include <utility/Exceptions.h>

using namespace ausaxs;

std::unique_ptr<hydrate::CullingStrategy> hydrate::factory::construct_culling_strategy(observer_ptr<data::Molecule> molecule) {
    return hydrate::factory::construct_culling_strategy(molecule, settings::hydrate::culling_strategy);
}

std::unique_ptr<hydrate::CullingStrategy> hydrate::factory::construct_culling_strategy(observer_ptr<data::Molecule> molecule, const settings::hydrate::CullingStrategy& choice) {
    switch (choice) {
        case settings::hydrate::CullingStrategy::CounterStrategy: 
            return std::make_unique<hydrate::CounterCulling>(molecule);
        case settings::hydrate::CullingStrategy::OutlierStrategy: 
            return std::make_unique<hydrate::OutlierCulling>(molecule);
        case settings::hydrate::CullingStrategy::NoStrategy:
            return std::make_unique<hydrate::NoCulling>(molecule);
        case settings::hydrate::CullingStrategy::RandomCounterStrategy:
            return std::make_unique<hydrate::RandomCulling<CounterCulling>>(molecule);
        case settings::hydrate::CullingStrategy::RandomOutlierStrategy:
            return std::make_unique<hydrate::RandomCulling<OutlierCulling>>(molecule);
        case settings::hydrate::CullingStrategy::RandomBodyCounterStrategy:
            return std::make_unique<hydrate::RandomCulling<BodyCounterCulling>>(molecule);
        default: 
            throw except::unknown_argument("grid::factory::construct_culling_strategy: Unkown CullingStrategy. Did you forget to add it to the switch statement?");
    }
}