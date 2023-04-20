#pragma once

#include <hydrate/culling/CounterCulling.h>
#include <hydrate/Grid.h>

#include <random>

namespace grid {
    /**
     * @brief Iterates through all placed water molecules, rejecting all but the nth, where n is determined from the desired number of water molecules. 
     */
    class RandomCulling : public CounterCulling {
        public:
            using CounterCulling::CounterCulling;
            ~RandomCulling() override = default;

            // runs in O(n) where n is the number of water molecules
            std::vector<Water> cull(std::vector<GridMember<Water>>& placed_water) const override {
                std::shuffle(placed_water.begin(), placed_water.end(), std::mt19937{std::random_device{}()}); // shuffle the molecules
                return CounterCulling::cull(placed_water);
            }
    };       
}