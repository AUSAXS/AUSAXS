// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hydrate/culling/CullingStrategy.h>

namespace ausaxs::hydrate {
    /**
     * @brief No culling will be performed.  
     */
    class NoCulling : public CullingStrategy {
        public:
            using CullingStrategy::CullingStrategy;
            ~NoCulling() override = default;

            void cull(std::span<grid::GridMember<data::Water>>& placed_water) const override;
    };       
}