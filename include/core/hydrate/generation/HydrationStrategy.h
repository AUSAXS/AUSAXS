#pragma once

#include <hydrate/Hydration.h>
#include <data/DataFwd.h>

namespace ausaxs::hydrate {
    class HydrationStrategy {
        public:
            HydrationStrategy() = default;
            virtual ~HydrationStrategy() = default;

            /**
             * @brief Generate a new hydration layer using the chosen strategy. 
             */
            virtual void hydrate() = 0;

            /**
             * @brief Determine if this strategy supports Body-specific hydration layers, or only a single global hydration layer.
             */
            virtual bool global() const = 0;
    };
}