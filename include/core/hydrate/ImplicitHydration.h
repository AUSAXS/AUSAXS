// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hydrate/Hydration.h>

namespace ausaxs::hydrate {
    class ImplicitHydration final : public Hydration {
        public:
            ImplicitHydration();

            ~ImplicitHydration();

            void clear() override;

            std::unique_ptr<Hydration> clone() const override;
    };
}