#pragma once

#include <data/DataFwd.h>
#include <hydrate/GridFwd.h>
#include <hydrate/detail/GridInternalFwd.h>

#include <vector>

namespace grid {    
    /**
     * @brief This class defines the strategy used to remove some of the water molecules. See its subclasses for more information on how this is done. 
     */
    class CullingStrategy {
        public:
            /**
             * @brief Constructor.
             * @param grid The Grid object to apply this Strategy to.
             */
            CullingStrategy(Grid* grid);

            /**
             * @brief Destructor.
             */
            virtual ~CullingStrategy();

            /**
             * @brief Cull the water molecules.
             * @return The remaining molecules after the culling.
             */
            virtual std::vector<data::record::Water> cull(std::vector<GridMember<data::record::Water>>& placed_water) const = 0;

            /**
             * @brief Set the desired number of water molecules after the culling. 
             * @param target_count The target number of water molecules. 
             */
            void set_target_count(unsigned int target_count);

        protected: 
            unsigned int target_count = 0; // The desired number of molecules after the culling.
            Grid* grid; // A reference to the grid used in Grid.
    };
}