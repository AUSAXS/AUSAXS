#pragma once

namespace constants {
    // Atom enum for a consistent way to denote the type of an atom. Example values: C, H, O, dummy
    enum class atom_t {
        H, He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, P, S, Cl, Ar, K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn, W, 
        M,      // This fake element is for compatibility with the GROMACS tip4p water model.
        dummy,  // Fake element for dummy atoms with variable properties. 
        unknown
    };

    // Atomic group enum for a consistent way to denote the type of an atomic group. Example values: CH, CH2, CH3, NH, NH2, OH, SH
    enum class atomic_group_t {
        // Note that any added group must also be added to the constants::symbols::get_atomic_group function to work with form factors. 
        CH, CH2, CH3, OH, NH, NH2, NH3, SH,
        unknown
    };
}