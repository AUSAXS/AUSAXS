#pragma once

#include <hydrate/HydrationStrategy.h>
#include <hydrate/grid/GridFwd.h>

#include <memory>

namespace hydrate {
    class GridBasedHydration : public HydrationStrategy {
        protected:
            virtual void hydrate_grid() = 0;
            std::unique_ptr<grid::Grid> grid;
    };
}