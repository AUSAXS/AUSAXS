// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

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
            bool expanded_across_symmetry = false; // Whether `waters` already covers the full symmetry-expanded assembly

            void clear() override;

            std::unique_ptr<Hydration> clone() const final override;
    };
}