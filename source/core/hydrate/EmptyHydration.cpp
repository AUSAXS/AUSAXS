// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hydrate/EmptyHydration.h>

using namespace ausaxs;

hydrate::EmptyHydration::EmptyHydration() = default;

hydrate::EmptyHydration::~EmptyHydration() = default;

void hydrate::EmptyHydration::clear() {}

std::unique_ptr<hydrate::Hydration> hydrate::EmptyHydration::clone() const {return std::make_unique<EmptyHydration>();}