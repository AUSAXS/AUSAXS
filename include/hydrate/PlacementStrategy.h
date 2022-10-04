#pragma once

class Grid;

#include <data/Water.h>
#include <hydrate/GridMember.h>

namespace grid {
    /**
     * @brief This class defines the strategy used to place water molecules. See its subclasses for more information on how this is done. 
     */
    class PlacementStrategy {
        public:
            /**
             * @brief Constructor. 
             * @param grid The Grid object to apply this Strategy to.
             */
            PlacementStrategy(Grid* grid) {this->grid = grid;}

            /**
             * @brief Destructor.
             */
            virtual ~PlacementStrategy() {}

            /**
             * @brief Place water molecules in the grid wherever possible.
             * @return A list of (binx, biny, binz) coordinates where the water molecules were placed.
             */
            virtual std::vector<GridMember<Water>> place() const = 0;

        protected: 
            Grid* grid; // A reference to the grid used in Grid.
    };
}