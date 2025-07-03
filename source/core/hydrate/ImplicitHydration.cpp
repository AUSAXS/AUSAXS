// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hydrate/ImplicitHydration.h>

using namespace ausaxs;

hydrate::ImplicitHydration::ImplicitHydration() = default;

hydrate::ImplicitHydration::~ImplicitHydration() = default;

void hydrate::ImplicitHydration::clear() {throw std::runtime_error("ImplicitHydration::clear: Not implemented.");}

std::unique_ptr<hydrate::Hydration> hydrate::ImplicitHydration::clone() const {throw std::runtime_error("ImplicitHydration::clone: Not implemented.");}