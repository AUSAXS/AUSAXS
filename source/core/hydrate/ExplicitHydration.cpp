// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hydrate/ExplicitHydration.h>
#include <data/atoms/Water.h>

using namespace ausaxs;

hydrate::ExplicitHydration::ExplicitHydration() = default;
hydrate::ExplicitHydration::ExplicitHydration(const std::vector<data::Water>& waters) : waters(waters) {}
hydrate::ExplicitHydration::ExplicitHydration(std::vector<data::Water>&& waters) : waters(std::move(waters)) {}
hydrate::ExplicitHydration::~ExplicitHydration() = default;

void hydrate::ExplicitHydration::clear() {
    waters.clear();
}

std::unique_ptr<hydrate::Hydration> hydrate::ExplicitHydration::clone() const {
    return std::make_unique<hydrate::ExplicitHydration>(*this);
}