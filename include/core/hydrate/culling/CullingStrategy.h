#pragma once

#include <data/DataFwd.h>
#include <utility/observer_ptr.h>
#include <grid/GridFwd.h>
#include <grid/detail/GridInternalFwd.h>

#include <vector>

namespace ausaxs::hydrate {    
    /**
     * @brief This class defines the strategy used to remove some of the water molecules from the grid. 
     *        See its subclasses for more information on how this is done. 
     */
    class CullingStrategy {
        public:
            CullingStrategy(observer_ptr<data::Molecule> grid);
            virtual ~CullingStrategy();

            /**
             * @brief Cull the water molecules.
             * @return A copy of the remaining molecules after the culling. They will also be stored in the grid itself. 
             */
            virtual std::vector<data::Water>& cull(std::vector<grid::GridMember<data::Water>>& placed_water) const = 0;

            /**
             * @brief Set the desired number of water molecules after the culling. 
             * @param target_count The target number of water molecules. 
             */
            void set_target_count(unsigned int target_count);

        protected: 
            unsigned int target_count = 0; // The desired number of molecules after the culling.
            observer_ptr<data::Molecule> molecule;
    };
}