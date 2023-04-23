#include <hydrate/culling/CullingFactory.h>
#include <hydrate/culling/CounterCulling.h>
#include <hydrate/culling/OutlierCulling.h>
#include <hydrate/culling/RandomCulling.h>
#include <utility/Exceptions.h>

std::unique_ptr<grid::CullingStrategy> grid::factory::construct_culling_strategy(Grid* grid, settings::grid::CullingStrategy choice) {
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