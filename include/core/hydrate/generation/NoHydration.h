#pragma once

#include <hydrate/generation/HydrationStrategy.h>

namespace ausaxs::hydrate {
    /**
     * @brief This strategy will not place any water molecules. 
     */
    class NoHydration : public HydrationStrategy {
        public:
            using HydrationStrategy::HydrationStrategy;
            ~NoHydration() override = default;

            bool global() const override {return false;}

            void hydrate() override;
        };
}