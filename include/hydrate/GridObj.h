#pragma once

#include <vector>
#include <math/Vector3.h>

/**
 * @brief A simple class that represents a 3D grid.
 *        Designed to make access more consistent.
 */
class GridObj {
    using T = char;

    public: 
        // The different states a grid cell can be in. 
        // Public enum to make it globally available.
        enum State {
            EMPTY = 0,
            H_CENTER = 'H',
            A_CENTER = 'A',
            H_AREA = 'h',
            A_AREA = 'a'
        };

        GridObj() {}

        GridObj(unsigned int x, unsigned int y, unsigned int z) : xdim(x), ydim(y), zdim(z), grid(x, std::vector<std::vector<T>>(y, std::vector<T>(z, EMPTY))) {}

        T& index(unsigned int x, unsigned int y, unsigned int z) {return grid[x][y][z];}
        const T& index(unsigned int x, unsigned int y, unsigned int z) const {return grid[x][y][z];}

        T& index(const Vector3<int>& v) {return index(v.x(), v.y(), v.z());}
        const T& index(const Vector3<int>& v) const {return index(v.x(), v.y(), v.z());}

        unsigned int xdim, ydim, zdim;
    private:
        std::vector<std::vector<std::vector<T>>> grid; // The actual grid. Datatype is char since we need at least four different values.
};