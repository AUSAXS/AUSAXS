#pragma once

#include <hydrate/Hydration.h>

namespace ausaxs::hydrate {
    /**
     * @brief Not implemented.
     */
    class ImplicitHydration : public Hydration {
        public:
            ImplicitHydration() = default;
            ~ImplicitHydration() = default;

            void clear() override {};
    };
}