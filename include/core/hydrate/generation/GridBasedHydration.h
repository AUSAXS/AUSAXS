#pragma once

#include <hydrate/generation/HydrationStrategy.h>
#include <hydrate/culling/CullingStrategy.h>
#include <utility/observer_ptr.h>
#include <grid/detail/GridInternalFwd.h>
#include <grid/GridFwd.h>

#include <span>
#include <memory>

namespace ausaxs::hydrate {
    class GridBasedHydration : public HydrationStrategy {
        public:
            GridBasedHydration(observer_ptr<data::Molecule> protein);
            GridBasedHydration(observer_ptr<data::Molecule> protein, std::unique_ptr<CullingStrategy> culling_strategy);
            virtual ~GridBasedHydration();

            void hydrate() override;

            void set_culling_strategy(std::unique_ptr<CullingStrategy> culling_strategy);

        protected:
            observer_ptr<data::Molecule> protein;

            virtual std::span<grid::GridMember<data::Water>> generate_explicit_hydration(std::span<grid::GridMember<data::AtomFF>> atoms) = 0;

            virtual void initialize();

        private:
            std::unique_ptr<CullingStrategy> culling_strategy;
    };
}