#pragma once

// forward declaration
class Grid;

// includes
#include "data/Hetatom.h"
#include "hydrate/Grid.h"

using std::vector, std::string, std::shared_ptr, std::unique_ptr;
using namespace ROOT;

/**
 * @brief This class defines the strategy used to remove some of the water molecules. See its subclasses for more information on how this is done. 
 */
class CullingStrategy {
public:
    /**
     * @brief Construct a new Culling Strategy object.
     * @param grid the Grid object to apply this Strategy to.
     */
    CullingStrategy(Grid* grid) : grid(grid) {}

    virtual ~CullingStrategy() {}

    /**
     * @brief Cull the water molecules.
     * @return The remaining molecules after the culling.
     */
    virtual vector<shared_ptr<Hetatom>> cull(vector<shared_ptr<Hetatom>> placed_water) const = 0;

    /**
     * @brief Set the desired number of molecules after the culling. 
     */
    void set_target_count(int target_count) {this->target_count = target_count;}

protected: 
    int target_count = 0; // The desired number of molecules after the culling.
    Grid* grid; // A reference to the grid used in Grid.
};