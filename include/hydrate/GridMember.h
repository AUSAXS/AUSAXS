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
        GridMember(T atom, std::vector<unsigned int> loc) : atom(atom), loc(loc) {}

        T atom; // the atom itself
        std::vector<unsigned int> loc; // the bin location of the Atom key
        bool expanded_volume = false; // whether the volume of this location has been expanded

        /**
         * @brief Equality operator. 
         *        Check if this GridMember atom is the same as the one provided. 
         * @param rhs External atom to compare against. 
         */
        bool operator==(const T& rhs) const {return atom == rhs;}

        /**
         * @brief Equality operator. 
         *        Check if this GridMember the same as another. 
         * @param rhs GridMember to compare against. 
         */
        bool operator==(const GridMember<T>& rhs) const {
            if (atom != rhs.atom) {return false;}
            if (loc != rhs.loc) {return false;}
            if (expanded_volume != rhs.expanded_volume) {return false;}
            return true;
        } 

        /**
         * @brief Assignment operator. 
         *        Assign new a new GridMember to this object. 
         * @param rhs The assigned GridMember. 
         */
        GridMember& operator=(const GridMember<T>& rhs) {
            atom = rhs.atom;
            loc = rhs.loc;
            expanded_volume = rhs.expanded_volume;
            return *this;
        }
    };
}