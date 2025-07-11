// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hydrate/culling/CullingStrategy.h>

#include <vector>

namespace ausaxs::hydrate {
    /**
     * @brief Iterates through all placed water molecules, rejecting all but the nth, where n is determined from the desired number of water molecules. 
     *        This is done on a per-body basis, allowing for more control over the hydration placement. 
     */
    class BodyCounterCulling : public CullingStrategy {
        public:
            using CullingStrategy::CullingStrategy;
            ~BodyCounterCulling() override = default;

            // runs in O(n) where n is the number of water molecules
            void cull(std::span<grid::GridMember<data::Water>>& placed_water) const override;

            void set_body_ratios(const std::vector<double>& body_ratios);

        private:
            std::vector<double> body_ratios;
    };       
}