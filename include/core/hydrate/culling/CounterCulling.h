// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include "grid/detail/GridMember.h"
#include <hydrate/culling/CullingStrategy.h>

namespace ausaxs::hydrate {
    /**
     * @brief Iterates through all placed water molecules, rejecting all but the nth, where n is determined from the desired number of water molecules. 
     */
    class CounterCulling : public CullingStrategy {
        public:
            using CullingStrategy::CullingStrategy;
            ~CounterCulling() override = default;

            // runs in O(n) where n is the number of water molecules
            void cull(std::span<grid::GridMember<data::Water>>& placed_water) const override;
    };       
}