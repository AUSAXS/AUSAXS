// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hydrate/Hydration.h>
#include <hydrate/ExplicitHydration.h>
#include <hydrate/ImplicitHydration.h>
#include <hydrate/EmptyHydration.h>

namespace ausaxs::hydrate {
    template<data::WaterVector T>
    std::unique_ptr<Hydration> Hydration::create(T&& data) {
        if constexpr (std::is_same_v<std::remove_cvref_t<T>, std::vector<data::Water>>) {
            return std::make_unique<ExplicitHydration>(std::forward<T>(data));
        } else {
            return std::make_unique<ImplicitHydration>(std::forward<T>(data));
        }
    }

    std::unique_ptr<Hydration> Hydration::create() {
        return std::make_unique<EmptyHydration>();
    }
}

template std::unique_ptr<ausaxs::hydrate::Hydration> ausaxs::hydrate::Hydration::create(std::vector<ausaxs::data::Water>&&);
template std::unique_ptr<ausaxs::hydrate::Hydration> ausaxs::hydrate::Hydration::create(std::vector<ausaxs::data::Water>&);
template std::unique_ptr<ausaxs::hydrate::Hydration> ausaxs::hydrate::Hydration::create(const std::vector<ausaxs::data::Water>&);