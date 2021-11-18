#pragma once

// forward declaration
class Grid;

// includes
#include <TVector3.h>
#include "data/Atom.h"
#include "data/Hetatom.h"
#include "Grid.h"

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

    /**
     * @brief Place water molecules in the grid wherever possible.
     * @param bounds the area to place molecules in. 
     * @return A list of (binx, biny, binz) coordinates where the water molecules were placed.
     */
    virtual vector<shared_ptr<Hetatom>> place(const vector<vector<int>> bounds) = 0;

protected: 
    Grid* grid; // A reference to the grid used in Grid.

private:
    /**
     * @brief Check if a water molecule can be placed at the given location. 
     * @param loc the location to be checked. 
     * @return True if this is an acceptable location, false otherwise.
     */
    virtual bool collision_check(const vector<int> loc) const = 0;
};