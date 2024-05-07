#pragma once

#include <hydrate/generation/HydrationStrategy.h>
#include <utility/observer_ptr.h>
#include <grid/GridFwd.h>
#include <settings/MoleculeSettings.h>

#include <memory>

namespace hydrate {
    namespace factory {
        std::unique_ptr<HydrationStrategy> construct_hydration_generator(observer_ptr<grid::Grid> grid);
        std::unique_ptr<HydrationStrategy> construct_hydration_generator(observer_ptr<grid::Grid> grid, const settings::hydrate::HydrationStrategy& choice);
    }
}