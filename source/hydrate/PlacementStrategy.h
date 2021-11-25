#pragma once

// forward declaration
class Grid;

// includes
#include <TVector3.h>
#include "data/Hetatom.h"
#include "hydrate/Grid.h"

using std::vector, std::string, std::shared_ptr, std::unique_ptr;
using namespace ROOT;

/**
 * @brief This class defines the strategy used to place water molecules. See its subclasses for more information on how this is done. 
 */
class PlacementStrategy {
public:
    /**
     * @brief Construct a new Placement Strategy object.
     * @param grid the Grid object to apply this Strategy to.
     */
    PlacementStrategy(Grid* grid) {this->grid = grid;}

    virtual ~PlacementStrategy() {}

    /**
     * @brief Place water molecules in the grid wherever possible.
     * @return A list of (binx, biny, binz) coordinates where the water molecules were placed.
     */
    virtual vector<shared_ptr<Hetatom>> place() const = 0;

protected: 
    Grid* grid; // A reference to the grid used in Grid.
};