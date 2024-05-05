/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <grid/culling/CullingFactory.h>
#include <grid/culling/CounterCulling.h>
#include <grid/culling/OutlierCulling.h>
#include <grid/culling/RandomCulling.h>
#include <settings/GridSettings.h>
#include <utility/Exceptions.h>

std::unique_ptr<grid::CullingStrategy> grid::factory::construct_culling_strategy(Grid* grid) {
    return grid::factory::construct_culling_strategy(grid, settings::grid::culling_strategy);
}

std::unique_ptr<grid::CullingStrategy> grid::factory::construct_culling_strategy(Grid* grid, const settings::grid::CullingStrategy& choice) {
    switch (choice) {
        case settings::grid::CullingStrategy::CounterStrategy: 
            return std::make_unique<CounterCulling>(grid);
        case settings::grid::CullingStrategy::OutlierStrategy: 
            return std::make_unique<OutlierCulling>(grid);
        case settings::grid::CullingStrategy::RandomStrategy:
            return std::make_unique<RandomCulling>(grid);
        default: 
            throw except::unknown_argument("grid::factory::construct_culling_strategy: Unkown CullingStrategy. Did you forget to add it to the switch statement?");
    }
}