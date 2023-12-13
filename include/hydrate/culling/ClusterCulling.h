#pragma once

#include <math/Vector3.h>
#include <hydrate/GridFwd.h>

#include <vector>

namespace grid {
    /**
     * @brief Iterates through all placed water molecules, rejecting all but the nth, where n is determined from the desired number of water molecules. 
     */
    class ClusterCulling {
        public:
            ClusterCulling(Grid* grid);
            ~ClusterCulling() = default;

            /**
             * @brief Remove atoms part of too small groups. 
             */
            std::vector<bool> remove_clusters(unsigned int min_group_size) const;

            /**
             * @brief Remove atoms with too few neighbours.
             *        Should be run before remove_clusters to avoid tendril-like groups surviving. 
             */
            std::vector<bool> remove_tendrils(unsigned int min_neighbours) const;

        private: 
            Grid* grid;
            std::vector<Vector3<int>> rot_bins_2ra;
            std::vector<Vector3<int>> rot_bins_3ra;

            void prepare_rotations();
    };       
}