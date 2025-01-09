#pragma once

#include <hydrate/Hydration.h>
#include <data/DataFwd.h>

#include <vector>

namespace ausaxs::hydrate {
    class ExplicitHydration : public Hydration {
        public:
            ExplicitHydration();
            ExplicitHydration(const std::vector<data::Water>& waters);
            ExplicitHydration(std::vector<data::Water>&& waters);
            ~ExplicitHydration() override;

            std::vector<data::Water> waters;

            void clear() override;

            std::unique_ptr<Hydration> clone() const final override;
    };
}