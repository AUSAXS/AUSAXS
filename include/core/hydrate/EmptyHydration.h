// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hydrate/Hydration.h>

namespace ausaxs::hydrate {
    class EmptyHydration final : public Hydration {
        public:
            EmptyHydration();

            ~EmptyHydration() override;

            void clear() override;

            std::unique_ptr<Hydration> clone() const override;
    };
}