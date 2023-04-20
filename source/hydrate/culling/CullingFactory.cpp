#include <hydrate/culling/CullingFactory.h>
#include <hydrate/culling/CounterCulling.h>
#include <hydrate/culling/OutlierCulling.h>
#include <hydrate/culling/RandomCulling.h>
#include <utility/Exceptions.h>

using namespace settings::grid;
using namespace settings::detail;

// extend the settings namespace with a new option
template<> std::string SmartOption<CullingStrategy>::get() const {return std::to_string(static_cast<int>(value));}
template<> void SmartOption<CullingStrategy>::set(const std::vector<std::string>& val) {
    value = static_cast<CullingStrategy>(std::stoi(val[0]));
}

SmartOption<CullingStrategy> culling_strategy(CullingStrategy::CounterStrategy);
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