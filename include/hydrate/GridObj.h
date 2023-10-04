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
            // The different states a grid cell can be in. 
            // Public enum to makes them globally available.
            enum State : char {
                EMPTY =    (1 << 0),
                W_CENTER = (1 << 1),    // Center of a hydrogen atom.
                A_CENTER = (1 << 2),    // Center of an atom.
                W_AREA =   (1 << 3),    // Area surrounding the center of a hydrogen atom.
                A_AREA =   (1 << 4),    // Area surrounding the center of an atom.
                VOLUME =   (1 << 5)     // Similar to EMPTY, but used to indicate that the cell is part of the volume of the protein. We need this to avoid double-counting cells when evaluating the volume.
            };

            /**
             * @brief Branchless function to check if a given bin is empty.
             */
            bool is_empty(unsigned int x, unsigned int y, unsigned int z) const {
                return index(x, y, z) & (EMPTY | VOLUME);
            }

            /**
             * @brief Branchless function to check if a given bin is part of an atomic volume.
             */
            bool is_atom_area(unsigned int x, unsigned int y, unsigned int z) const {
                return index(x, y, z) & (A_AREA | VOLUME);
            }

            /**
             * @brief Branchless function to check if a given bin is part of a water volume.
             */
            bool is_water_area(unsigned int x, unsigned int y, unsigned int z) const {
                return index(x, y, z) & W_AREA;
            }

            /**
             * @brief Branchless function to check if a given bin is the center of an atom.
             */
            bool is_atom_center(unsigned int x, unsigned int y, unsigned int z) const {
                return index(x, y, z) & A_CENTER;
            }

            /**
             * @brief Branchless function to check if a given bin is the center of a water molecule.
             */
            bool is_water_center(unsigned int x, unsigned int y, unsigned int z) const {
                return index(x, y, z) & W_CENTER;
            }

            GridObj() = default;

            GridObj(unsigned int x, unsigned int y, unsigned int z);

            State& index(unsigned int x, unsigned int y, unsigned int z);
            const State& index(unsigned int x, unsigned int y, unsigned int z) const;

            State& index(const Vector3<int>& v);
            const State& index(const Vector3<int>& v) const;

            unsigned int xdim, ydim, zdim;
        private:
            std::vector<std::vector<std::vector<State>>> grid; // The actual grid. Datatype is char since we need at least four different values.
    };
}