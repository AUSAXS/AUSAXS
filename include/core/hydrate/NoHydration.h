#pragma once

#include <hydrate/Hydration.h>

namespace ausaxs::hydrate {
    class NoHydration : public Hydration {
        public:
            NoHydration() = default;
            ~NoHydration() override = default;
            void clear() override {}
            std::unique_ptr<Hydration> clone() const final override {return std::make_unique<NoHydration>();}
    };
}