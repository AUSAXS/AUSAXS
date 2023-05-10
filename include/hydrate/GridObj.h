#pragma once

#include <utility/Concepts.h>

#include <vector>

template<numeric T> class Vector3;
namespace grid {
    /**
     * @brief A simple class that represents a 3D grid.
     *        Designed to make access more consistent.
     */
    class GridObj {
        public: 
            using DATATYPE = char;

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

            GridObj(unsigned int x, unsigned int y, unsigned int z);

            DATATYPE& index(unsigned int x, unsigned int y, unsigned int z);
            const DATATYPE& index(unsigned int x, unsigned int y, unsigned int z) const;

            DATATYPE& index(const Vector3<int>& v);
            const DATATYPE& index(const Vector3<int>& v) const;

            unsigned int xdim, ydim, zdim;
        private:
            std::vector<std::vector<std::vector<DATATYPE>>> grid; // The actual grid. Datatype is char since we need at least four different values.
    };
}