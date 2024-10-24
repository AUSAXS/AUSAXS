#pragma once

#include <hydrate/Hydration.h>
#include <data/DataFwd.h>

#include <vector>

namespace ausaxs::hydrate {
    class ExplicitHydration : public Hydration {
        public:
            ExplicitHydration();
            ExplicitHydration(const std::vector<data::record::Water>& waters);
            ExplicitHydration(std::vector<data::record::Water>&& waters);
            ~ExplicitHydration() override;

            std::vector<data::record::Water> waters;

            void clear() override;
    };
}