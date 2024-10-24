#pragma once

#include <hydrate/Hydration.h>
#include <data/DataFwd.h>

#include <memory>

namespace ausaxs::hydrate {
    class HydrationStrategy {
        public:
            HydrationStrategy() = default;
            virtual ~HydrationStrategy() = default;

            virtual std::unique_ptr<Hydration> hydrate() = 0;
    };
}