#pragma once

#include <hydrate/generation/GridBasedHydration.h>
#include <data/DataFwd.h>
#include <math/MathFwd.h>

namespace hydrate {
    /**
     * @brief This strategy mimics the model proposed by PepsiSAXS 10.1107/S2059798317005745 for generating explicit hydration sites.
     */
    class PepsiHydration : public GridBasedHydration {
        public:
            PepsiHydration(observer_ptr<data::Molecule> protein);
            ~PepsiHydration() override;

            std::vector<data::record::Water> generate_explicit_hydration() override;

        private:
            void modified_expand_volume(grid::GridMember<data::record::Atom>& atom);
            void initialize() override;
    };
}