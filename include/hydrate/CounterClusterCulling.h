#pragma once

#include <hydrate/ClusterCullingStrategy.h>
#include <hydrate/Grid.h>

namespace grid {
    /**
     * @brief Iterates through all placed water molecules, rejecting all but the nth, where n is determined from the desired number of water molecules. 
     */
    class CounterClusterCulling : public ClusterCullingStrategy {
        public:
            CounterClusterCulling(Grid* grid);
            ~CounterClusterCulling() override {}

            // runs in O(n) where n is the number of water molecules
            std::vector<bool> cull(unsigned int min_group_size) const override;
        
        private: 
            std::vector<Vector3<int>> rot_bins_2ra;
            std::vector<Vector3<int>> rot_bins_3ra;

            void prepare_rotations();
    };       
}