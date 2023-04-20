#include <hydrate/placement/PlacementFactory.h>
#include <hydrate/placement/JanPlacement.h>
#include <hydrate/placement/RadialPlacement.h>
#include <hydrate/placement/AxesPlacement.h>
#include <utility/Exceptions.h>

using namespace settings::grid;
using namespace settings::detail;

// extend the settings namespace with a new option
template<> std::string SmartOption<PlacementStrategy>::get() const {return std::to_string(static_cast<int>(value));}
template<> void SmartOption<PlacementStrategy>::set(const std::vector<std::string>& val) {
    value = static_cast<PlacementStrategy>(std::stoi(val[0]));
}

SmartOption<PlacementStrategy> placement_strategy(PlacementStrategy::RadialStrategy);
std::unique_ptr<grid::PlacementStrategy> grid::factory::construct_placement_strategy(Grid* grid, settings::grid::PlacementStrategy choice) {
    switch (choice) {
        case settings::grid::PlacementStrategy::AxesStrategy: 
            return std::make_unique<AxesPlacement>(grid);
        case settings::grid::PlacementStrategy::RadialStrategy:
            return std::make_unique<RadialPlacement>(grid);
        case settings::grid::PlacementStrategy::JanStrategy: 
            return std::make_unique<JanPlacement>(grid);
        default: 
            throw except::unknown_argument("grid::factory::construct_placement_strategy: Unkown PlacementStrategy. Did you forget to add it to the switch statement?");
    }
}
