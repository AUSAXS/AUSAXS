#pragma once

#include <hydrate/generation/HydrationStrategy.h>
#include <hydrate/culling/CullingStrategy.h>
#include <utility/observer_ptr.h>
#include <grid/GridFwd.h>

namespace hydrate {
    class GridBasedHydration : public HydrationStrategy {
        public:
            GridBasedHydration(observer_ptr<grid::Grid> grid);
            virtual ~GridBasedHydration();

            std::unique_ptr<Hydration> hydrate() override;

        protected:
            observer_ptr<grid::Grid> grid;

            virtual std::vector<data::record::Water> generate_explicit_hydration() = 0;

        private:
            std::unique_ptr<CullingStrategy> culling_strategy;
    };
}