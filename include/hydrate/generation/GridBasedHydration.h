#pragma once

#include <hydrate/generation/HydrationStrategy.h>
#include <hydrate/culling/CullingStrategy.h>
#include <utility/observer_ptr.h>
#include <grid/GridFwd.h>

namespace hydrate {
    class GridBasedHydration : public HydrationStrategy {
        public:
            GridBasedHydration(observer_ptr<data::Molecule> protein);
            virtual ~GridBasedHydration();

            std::unique_ptr<Hydration> hydrate() override;

        protected:
            observer_ptr<data::Molecule> protein;

            virtual std::vector<data::record::Water> generate_explicit_hydration() = 0;

            virtual void initialize();

        private:
            std::unique_ptr<CullingStrategy> culling_strategy;
    };
}