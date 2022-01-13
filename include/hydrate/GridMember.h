#pragma once

#include <vector>

template<typename T>
struct GridMember {
    GridMember() {}
    GridMember(const GridMember<T>& gm) : atom(gm.atom), loc(gm.loc), expanded_volume(gm.expanded_volume) {}
    GridMember(const GridMember<T>&& gm) noexcept : atom(std::move(gm.atom)), loc(std::move(gm.loc)), expanded_volume(std::move(gm.expanded_volume)) {}
    GridMember(T atom, std::vector<int> loc) : atom(atom), loc(loc) {}

    T atom; // the atom itself
    std::vector<int> loc; // the bin location of the Atom key
    bool expanded_volume = false; // whether the volume of this location has been expanded

    // the two operator overloads makes this struct act just like a vector, but with an additional bool available when needed. 
    int operator[](int index) {return loc[index];}
    int operator[](int index) const {return loc[index];}

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
};