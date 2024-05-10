#pragma once

#include <hydrate/generation/HydrationStrategy.h>

namespace hydrate {
    /**
     * @brief This strategy will not place any water molecules. 
     */
    class NoHydration : public HydrationStrategy {
        public:
            using HydrationStrategy::HydrationStrategy;
            ~NoHydration() override = default;

            std::unique_ptr<Hydration> hydrate() override;
        };
}