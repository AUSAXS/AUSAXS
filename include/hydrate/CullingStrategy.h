#pragma once

// forward declaration
class Grid;

// includes
#include "data/Hetatom.h"
#include "hydrate/Grid.h"
#include "hydrate/GridMember.h"

using std::vector, std::string, std::shared_ptr, std::unique_ptr;

/**
 * @brief This class defines the strategy used to remove some of the water molecules. See its subclasses for more information on how this is done. 
 */
class CullingStrategy {
public:
    /**
     * @brief Constructor.
     * @param grid The Grid object to apply this Strategy to.
     */
    CullingStrategy(Grid* grid) : grid(grid) {}

    /**
     * @brief Destructor.
     */
    virtual ~CullingStrategy() = default;

    /**
     * @brief Cull the water molecules.
     * @return The remaining molecules after the culling.
     */
    virtual vector<Hetatom> cull(vector<GridMember<Hetatom>>& placed_water) const = 0;

    /**
     * @brief Set the desired number of water molecules after the culling. 
     * @param target_count The target number of water molecules. 
     */
    void set_target_count(size_t target_count) {this->target_count = target_count;}

protected: 
    size_t target_count = 0; // The desired number of molecules after the culling.
    Grid* grid; // A reference to the grid used in Grid.
};