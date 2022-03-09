#pragma once

#include <em/ImageStack.h>

using std::vector, std::list, std::string;

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
         * @brief Cull the water molecules.
         * @return The remaining molecules after the culling.
         */
        virtual vector<Atom> cull(list<Atom>& atoms) const = 0;

        /**
         * @brief Set the desired number of water molecules after the culling. 
         * @param target_count The target number of water molecules. 
         */
        void set_target_count(size_t target_count) {this->target_count = target_count;}

    protected: 
        size_t target_count = 0; // The desired number of molecules after the culling.
        Grid* grid; // A reference to the grid used in Grid.
    };
}