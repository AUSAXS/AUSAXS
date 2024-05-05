#pragma once

#include <hydrate/GridBasedHydration.h>

namespace hydrate {
    class RadialHydration : public GridBasedHydration {
        public:
            RadialHydration() = default;
            virtual ~RadialHydration() = default;

            void hydrate() override;

        protected:
            void hydrate_grid() override;
    };
}