#pragma once

#include <utility/Concepts.h>
#include <math/MathFwd.h>
#include <container/Container3D.h>

#include <vector>

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
             * @brief Branchless function to check if a given bin is part of a volume. This means the bin is VOLUME.
             */
            bool is_volume(unsigned int x, unsigned int y, unsigned int z) const {return index(x, y, z) & VOLUME;}

            /**
             * @brief Branchless function to check if a given bin is part of a volume. This means the bin is VOLUME.
             */
            bool is_volume(State s) const {return s & VOLUME;}

            /**
             * @brief Branchless function to check if a given bin is empty or part of a volume. This means the bin is EMPTY or VOLUME.
             */
            bool is_empty_or_volume(unsigned int x, unsigned int y, unsigned int z) const {return index(x, y, z) & (EMPTY | VOLUME);}

            /**
             * @brief Branchless function to check if a given bin is empty or part of a volume. This means the bin is EMPTY or VOLUME.
             */
            bool is_empty_or_volume(State s) const {return s & (EMPTY | VOLUME);}

            /**
             * @brief Branchless function to check if a given bin is empty. Empty means the bin is EMPTY. 
             */
            bool is_empty(unsigned int x, unsigned int y, unsigned int z) const {return index(x, y, z) & EMPTY;}

            /**
             * @brief Branchless function to check if a given bin is empty. Empty means the bin is EMPTY. 
             */
            bool is_empty(State s) const {return s & EMPTY;}

            /**
             * @brief Branchless function to check if a given bin is part of an atomic volume. This means the bin is A_AREA. 
             */
            bool is_atom_area(unsigned int x, unsigned int y, unsigned int z) const {return index(x, y, z) & A_AREA;}

            /**
             * @brief Branchless function to check if a given bin is part of an atomic volume. This means the bin is A_AREA. 
             */
            bool is_atom_area(State s) const {return s & A_AREA;}

            /**
             * @brief Branchless function to check if a given bin is part of an atomic volume. This means the bin is either A_AREA or VOLUME.
             */
            bool is_atom_area_or_volume(unsigned int x, unsigned int y, unsigned int z) const {return index(x, y, z) & (A_AREA | VOLUME);}

            /**
             * @brief Branchless function to check if a given bin is part of an atomic volume. This means the bin is either A_AREA or VOLUME.
             */
            bool is_atom_area_or_volume(State s) const {return s & (A_AREA | VOLUME);}

            /**
             * @brief Branchless function to check if a given bin is part of a water volume. This means the bin is W_AREA.
             */
            bool is_water_area(unsigned int x, unsigned int y, unsigned int z) const {return index(x, y, z) & W_AREA;}

            /**
             * @brief Branchless function to check if a given bin is part of a water volume. This means the bin is W_AREA.
             */
            bool is_water_area(State s) const {return s & W_AREA;}

            /**
             * @brief Branchless function to check if a given bin is the center of an atom. This means the bin is A_CENTER.
             */
            bool is_atom_center(unsigned int x, unsigned int y, unsigned int z) const {return index(x, y, z) & A_CENTER;}

            /**
             * @brief Branchless function to check if a given bin is the center of an atom. This means the bin is A_CENTER.
             */
            bool is_atom_center(State s) const {return s & A_CENTER;}

            /**
             * @brief Branchless function to check if a given bin is the center of a water molecule. This means the bin is W_CENTER. 
             */
            bool is_water_center(unsigned int x, unsigned int y, unsigned int z) const {return index(x, y, z) & W_CENTER;}

            /**
             * @brief Branchless function to check if a given bin is the center of a water molecule. This means the bin is W_CENTER. 
             */
            bool is_water_center(State s) const {return s & W_CENTER;}

            GridObj() = default;

            GridObj(unsigned int x, unsigned int y, unsigned int z);

            State& index(unsigned int x, unsigned int y, unsigned int z);
            const State& index(unsigned int x, unsigned int y, unsigned int z) const;

            State& index(const Vector3<int>& v);
            const State& index(const Vector3<int>& v) const;

            unsigned int xdim, ydim, zdim;
        private:
            container::Container3D<State> grid; // The actual grid. Datatype is char since we need at least four different values.
    };
}