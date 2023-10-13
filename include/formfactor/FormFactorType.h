#pragma once

namespace form_factor {
    // The form factor type of an atom. This is intended to be used as an index for best performance.
    enum class form_factor_t {
        C,                  // neutral carbon
        CH,                 // neutral carbon with hydrogen
        CH2,                // neutral carbon with two hydrogens
        CH3,                // neutral carbon with three hydrogens
        N,                  // neutral nitrogen
        NH,                 // neutral nitrogen with hydrogen
        NH2,                // neutral nitrogen with two hydrogens
        NH3,                // neutral nitrogen with three hydrogens
        O,                  // neutral oxygen
        OH,                 // neutral oxygen with hydrogen
        S,                  // neutral sulfur
        SH,                 // neutral sulfur with hydrogen
        OTHER,              // all other atoms
        EXCLUDED_VOLUME,    // excluded volume
        COUNT,              // this will have the numerical value of the number of form factor types, and can thus be used to allocate arrays
    };
}