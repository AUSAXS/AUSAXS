/*
This software is distributed under the GNU General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <hydrate/placement/PlacementFactory.h>
#include <hydrate/placement/JanPlacement.h>
#include <hydrate/placement/RadialPlacement.h>
#include <hydrate/placement/AxesPlacement.h>
#include <settings/GridSettings.h>
#include <utility/Exceptions.h>
#include <math/Vector3.h>

std::unique_ptr<grid::PlacementStrategy> grid::factory::construct_placement_strategy(Grid* grid) {
    return grid::factory::construct_placement_strategy(grid, settings::grid::placement_strategy);
}

std::unique_ptr<grid::PlacementStrategy> grid::factory::construct_placement_strategy(Grid* grid, const settings::grid::PlacementStrategy& choice) {
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
