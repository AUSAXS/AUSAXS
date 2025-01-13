#pragma once

#include <utility/Concepts.h>
#include <math/MathFwd.h>
#include <container/Container3D.h>

namespace ausaxs::grid::detail {
    /**
        * @brief A simple enum to represent the different states of a grid cell. 
        *        We use bit flags to allow more efficient bitwise operations on them. 
        */
    enum State : uint8_t {
        EMPTY       = (1 << 0), // empty bin
        VOLUME      = (1 << 1), // bin is part of atomic volume outside its vdw radius
        VACUUM      = (1 << 2), 
        A_CENTER    = (1 << 3), // bin is the center of an atom
        A_AREA      = (1 << 4), // bin is part of an atomic volume inside its vdw radius
        W_CENTER    = (1 << 5), // bin is the center of a water molecule
        W_AREA      = (1 << 6), // bin is part of a water volume inside its vdw radius
        RESERVED_1  = (1 << 7),
    };

    /**
        * @brief A simple class that represents a 3D grid.
        *        Designed to make access more consistent.
        */
    class GridObj : public container::Container3D<detail::State> {
        public: 
            GridObj() = default;

            GridObj(unsigned int x, unsigned int y, unsigned int z);

            /**
                * @brief Branchless function to check if a given bin is part of a volume. This means the bin is VOLUME.
                */
            bool is_volume(State s) const;

            // @copydoc is_volume(State s) const
            bool is_volume(unsigned int x, unsigned int y, unsigned int z) const;

            /**
                * @brief Branchless function to check if a given bin is empty or part of a volume. This means the bin is EMPTY or VOLUME.
                */
            bool is_empty_or_volume(State s) const;

            // @copydoc is_empty_or_volume(State) const
            bool is_empty_or_volume(unsigned int x, unsigned int y, unsigned int z) const;

            /**
                * @brief Branchless function to check if a given bin is empty or part of the hydration shell. 
                *        This means the bin is EMPTY, W_AREA, or W_CENTER.
                */
            bool is_empty_or_water(State s) const;

            // @copydoc is_empty_or_water(State) const
            bool is_empty_or_water(unsigned int x, unsigned int y, unsigned int z) const;

            /**
                * @brief Branchless function to check if a given bin is empty, part of a volume, or part of the hydration shell. 
                *      This means the bin is EMPTY, VOLUME, W_AREA, or W_CENTER.
                */
            bool is_empty_or_volume_or_water(State s) const;

            // @copydoc is_empty_or_volume_or_water(State) const
            bool is_empty_or_volume_or_water(unsigned int x, unsigned int y, unsigned int z) const;

            /**
                * @brief Branchless function to check if a given bin is empty. Empty means the bin is EMPTY. 
                */
            bool is_empty(State s) const;

            // @copydoc is_empty(State) const
            bool is_empty(unsigned int x, unsigned int y, unsigned int z) const;

            /**
                * @brief Branchless function to check if a given bin is part of an atomic volume. This means the bin is A_AREA. 
                */
            bool is_atom_area(State s) const;

            // @copydoc is_atom_area(State) const
            bool is_atom_area(unsigned int x, unsigned int y, unsigned int z) const;

            /**
                * @brief Branchless function to check if a given bin is part of an atomic volume. This means the bin is either A_AREA or VOLUME.
                */
            bool is_atom_area_or_volume(State s) const;

            // @copydoc is_atom_area_or_volume(State) const
            bool is_atom_area_or_volume(unsigned int x, unsigned int y, unsigned int z) const;

            /**
                * @brief Branchless function to check if a given bin is part of a water volume. This means the bin is W_AREA.
                */
            bool is_water_area(State s) const;

            // @copydoc is_water_area(State) const
            bool is_water_area(unsigned int x, unsigned int y, unsigned int z) const;

            /**
                * @brief Branchless function to check if a given bin is the center of an atom. This means the bin is A_CENTER.
                */
            bool is_atom_center(State s) const;

            // @copydoc is_atom_center(State) const;
            bool is_atom_center(unsigned int x, unsigned int y, unsigned int z) const;

            /**
                * @brief Branchless function to check if a given bin is the center of a water molecule. This means the bin is W_CENTER. 
                */
            bool is_water_center(State s) const;

            // @copydoc is_water_center(State) const;
            bool is_water_center(unsigned int x, unsigned int y, unsigned int z) const;

            using container::Container3D<State>::index;
            State& index(const Vector3<int>& v);
            const State& index(const Vector3<int>& v) const;
    };

    constexpr grid::detail::State operator|(grid::detail::State lhs, grid::detail::State rhs) {
        return static_cast<grid::detail::State>(
            static_cast<std::underlying_type<grid::detail::State>::type>(lhs) |
            static_cast<std::underlying_type<grid::detail::State>::type>(rhs)
        );
    }

    constexpr grid::detail::State operator&(grid::detail::State lhs, grid::detail::State rhs) {
        return static_cast<grid::detail::State>(
            static_cast<std::underlying_type<grid::detail::State>::type>(lhs) &
            static_cast<std::underlying_type<grid::detail::State>::type>(rhs)
        );
    }

    constexpr grid::detail::State operator|=(grid::detail::State& lhs, grid::detail::State rhs) {
        lhs = lhs | rhs;
        return lhs;
    }

    constexpr grid::detail::State operator&=(grid::detail::State& lhs, grid::detail::State rhs) {
        lhs = lhs & rhs;
        return lhs;
    }

    constexpr grid::detail::State operator~(grid::detail::State s) {
        return static_cast<grid::detail::State>(~static_cast<std::underlying_type<grid::detail::State>::type>(s));
    }
}