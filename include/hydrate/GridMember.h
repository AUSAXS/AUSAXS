#pragma once

#include <vector>

namespace grid {
    /**
     * @brief \class GridMember
     * 
     * This class is used internally in Grid for storing all information about a particular member atom. 
     */
    template<typename T>
    struct GridMember {
        /**
         * @brief Default constructor.
         */
        GridMember() {}

        /**
         * @brief Copy constructor.
         */
        GridMember(const GridMember<T>& gm) : atom(gm.atom), loc(gm.loc), expanded_volume(gm.expanded_volume) {}

        /**
         * @brief Move constructor.
         */
        GridMember(const GridMember<T>&& gm) noexcept : atom(std::move(gm.atom)), loc(std::move(gm.loc)), expanded_volume(std::move(gm.expanded_volume)) {}

        /**
         * @brief Constructor.
         * @param atom The atom itself. 
         * @param loc The grid location of the atom. 
         */
        GridMember(T atom, Vector3<int> loc) : atom(atom), loc(std::move(loc)) {}

        // The atom itself
        T atom;

        // The bin location of the Atom key. Although this should never be negative, we don't use unsigned ints because it's a mess with implicit conversion errors.
        Vector3<int> loc; 

        // Whether the volume of this location has been expanded or not.
        bool expanded_volume = false; 

        bool operator==(const T& rhs) const {return atom == rhs;}
        bool operator==(const GridMember<T>& rhs) const {
            if (atom != rhs.atom) {return false;}
            if (loc != rhs.loc) {return false;}
            if (expanded_volume != rhs.expanded_volume) {return false;}
            return true;
        } 

        GridMember& operator=(const GridMember<T>& rhs) {
            atom = rhs.atom;
            loc = rhs.loc;
            expanded_volume = rhs.expanded_volume;
            return *this;
        }

        GridMember& operator=(GridMember<T>&& rhs) {
            atom = std::move(rhs.atom);
            loc = std::move(rhs.loc);
            expanded_volume = std::move(rhs.expanded_volume);
            return *this;
        }
    };
}