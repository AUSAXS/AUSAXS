#include <grid/expansion/GridExpanderFactory.h>
#include <grid/expansion/MinimalExpander.h>
#include <grid/expansion/FullExpander.h>

#include <stdexcept>

using namespace ausaxs::grid::expander;

std::unique_ptr<GridExpander> factory::construct(observer_ptr<Grid> grid, settings::grid::exv::Expansion choice) {
    switch (choice) {
        case settings::grid::exv::Expansion::Minimal: 
            return std::make_unique<MinimalExpander>(grid);
        case settings::grid::exv::Expansion::Full:
            return std::make_unique<FullExpander>(grid);
        default:
            throw std::invalid_argument(
                "grid::expander::factory::construct: Unknown expansion choice. "
                "Did you forget to add it to the switch statement?"
            );
    }
}

std::unique_ptr<GridExpander> factory::construct(observer_ptr<Grid> grid) {
    return construct(grid, settings::grid::exv::expansion_strategy);
}