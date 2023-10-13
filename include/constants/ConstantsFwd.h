#pragma once

namespace constants {
    // Atom enum for a consistent way to denote the type of an atom
    enum class atom_t {
        e, H, He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, P, S, Cl, Ar, K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn, W, 
        M,      // This fake element is for compatibility with the GROMACS tip4p water model.
        dummy,  // Fake element for dummy atoms with variable properties. 
        unknown
    };
}