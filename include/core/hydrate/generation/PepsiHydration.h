// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hydrate/generation/GridBasedHydration.h>
#include <data/DataFwd.h>
#include <math/MathFwd.h>

namespace ausaxs::hydrate {
    /**
     * @brief This strategy mimics the model proposed by PepsiSAXS 10.1107/S2059798317005745 for generating explicit hydration sites.
     */
    class PepsiHydration : public GridBasedHydration {
        public:
            PepsiHydration(observer_ptr<data::Molecule> protein);
            PepsiHydration(observer_ptr<data::Molecule> protein, std::unique_ptr<CullingStrategy> culling_strategy);
            ~PepsiHydration() override;

            std::span<grid::GridMember<data::Water>> generate_explicit_hydration(std::span<grid::GridMember<data::AtomFF>> atoms) override;

            bool global() const override {return false;}

        private:
            void initialize() override;
    };
}