#pragma once

namespace hydrate {
    class HydrationStrategy {
        public:
            HydrationStrategy() = default;
            virtual ~HydrationStrategy() = default;

            virtual void hydrate() = 0;
    };
}