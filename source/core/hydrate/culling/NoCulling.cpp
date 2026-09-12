// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hydrate/culling/NoCulling.h>

#include <grid/detail/GridMember.h>

void ausaxs::hydrate::NoCulling::cull(std::span<grid::GridMember<data::Water>>& /*placed_water*/) const {}