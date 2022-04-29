#pragma once

#include <vector>
#include <list>

#include <data/Atom.h>
#include <settings.h>

namespace em {
    /**
     * @brief This class defines the strategy used to remove some of the water molecules. See its subclasses for more information on how this is done. 
     */
    class CullingStrategy {
        public:
            /**
             * @brief Constructor.
             */
            CullingStrategy() {}

            /**
             * @brief Destructor.
             */
            virtual ~CullingStrategy() = default;

            /**
             * @brief Cull the atoms.
             * 
             * @return The remaining atoms after the culling.
             */
            virtual std::vector<Atom> cull(std::list<Atom>& atoms) const = 0;

            /**
             * @brief Set the desired number of atoms after the culling. 
             * 
             * @param target_count The target number of atoms. 
             */
            void set_target_count(size_t target_count) {this->target_count = target_count;}

        protected: 
            // unsigned int target_count = setting::em::max_atoms; // The desired number of molecules after the culling.
            unsigned int target_count = 50000; // The desired number of molecules after the culling.
    };
}