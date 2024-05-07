/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <hydrate/culling/CullingFactory.h>
#include <hydrate/culling/CounterCulling.h>
#include <hydrate/culling/OutlierCulling.h>
#include <hydrate/culling/RandomCulling.h>
#include <settings/GridSettings.h>
#include <utility/Exceptions.h>

std::unique_ptr<hydrate::CullingStrategy> hydrate::factory::construct_culling_strategy(observer_ptr<grid::Grid> grid) {
    return hydrate::factory::construct_culling_strategy(grid, settings::hydrate::culling_strategy);
}

std::unique_ptr<hydrate::CullingStrategy> hydrate::factory::construct_culling_strategy(observer_ptr<grid::Grid> grid, const settings::hydrate::CullingStrategy& choice) {
    switch (choice) {
        case settings::hydrate::CullingStrategy::CounterStrategy: 
            return std::make_unique<hydrate::CounterCulling>(grid);
        case settings::hydrate::CullingStrategy::OutlierStrategy: 
            return std::make_unique<hydrate::OutlierCulling>(grid);
        case settings::hydrate::CullingStrategy::RandomStrategy:
            return std::make_unique<hydrate::RandomCulling>(grid);
        default: 
            throw except::unknown_argument("grid::factory::construct_culling_strategy: Unkown CullingStrategy. Did you forget to add it to the switch statement?");
    }
}