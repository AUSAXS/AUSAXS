#pragma once

class Grid;

#include <vector>

#include <data/Hetatom.h>
#include <hydrate/GridMember.h>

namespace grid {
    /**
     * @brief This class defines the strategy used to count & remove floating clusters of water molecules. See its subclasses for more information on how this is done.
     */
    class ClusterCullingStrategy {
        public:
            /**
             * @brief Constructor.
             * @param grid The Grid object to apply this Strategy to.
             */
            ClusterCullingStrategy(Grid* grid) : grid(grid) {}

            /**
             * @brief Destructor.
             */
            virtual ~ClusterCullingStrategy() = default;

            virtual std::vector<bool> cull(unsigned int min_group_size) const = 0;

        protected: 
            Grid* grid; // A reference to the grid used in Grid.
    };
}