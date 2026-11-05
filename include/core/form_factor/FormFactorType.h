// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

namespace ausaxs::form_factor {
    // The form factor type of an atom. This is intended to be used as an index for best performance.
    // NOLINTNEXTLINE(readability-enum-initial-value)
    enum class form_factor_t : int {
        EXCLUDED_VOLUME,    // excluded volume
        WATER,              // water
        OH = WATER,         // neutral oxygen with hydrogen. Must directly follow WATER, since the enumerators after it count on from its value

        // BIOMOLECULAR FORM FACTORS
        C,                  // neutral carbon
        CH,                 // neutral carbon with hydrogen
        CH2,                // neutral carbon with two hydrogens
        CH3,                // neutral carbon with three hydrogens
        N,                  // neutral nitrogen
        NH,                 // neutral nitrogen with hydrogen
        NH2,                // neutral nitrogen with two hydrogens
        NH3,                // neutral nitrogen with three hydrogens
        O,                  // neutral oxygen
        S,                  // neutral sulfur
        SH,                 // neutral sulfur with hydrogen

        // ATOMIC FORM FACTORS
        H,                  // neutral hydrogen
        D,                  // neutral deuterium
        He,                 // neutral helium
        Li,                 // neutral lithium
        Be,                 // neutral beryllium
        B,                  // neutral boron
        F,                  // neutral fluorine
        Ne,                 // neutral neon
        Na,                 // neutral sodium
        Mg,                 // neutral magnesium
        Al,                 // neutral aluminum
        Si,                 // neutral silicon
        P,                  // neutral phosphorus
        Cl,                 // neutral chlorine
        Ar,                 // neutral argon
        K,                  // neutral potassium
        Ca,                 // neutral calcium
        Sc,                 // neutral scandium
        Ti,                 // neutral titanium
        V,                  // neutral vanadium
        Cr,                 // neutral chromium
        Mn,                 // neutral manganese
        Fe,                 // neutral iron
        Co,                 // neutral cobalt
        Ni,                 // neutral nickel
        Cu,                 // neutral copper
        Zn,                 // neutral zinc
        Ga,                 // neutral gallium
        Ge,                 // neutral germanium
        As,                 // neutral arsenic
        Se,                 // neutral selenium
        Br,                 // neutral bromine
        Kr,                 // neutral krypton
        Rb,                 // neutral rubidium
        Sr,                 // neutral strontium
        Y,                  // neutral yttrium
        Zr,                 // neutral zirconium
        Nb,                 // neutral niobium
        Mo,                 // neutral molybdenum
        Tc,                 // neutral technetium
        Ru,                 // neutral ruthenium
        Rh,                 // neutral rhodium
        Pd,                 // neutral palladium
        Ag,                 // neutral silver
        Cd,                 // neutral cadmium
        In,                 // neutral indium
        Sn,                 // neutral tin
        Sb,                 // neutral antimony
        Te,                 // neutral tellurium
        I,                  // neutral iodine
        Xe,                 // neutral xenon
        Cs,                 // neutral cesium
        Ba,                 // neutral barium
        Hf,                 // neutral hafnium
        Ta,                 // neutral tantalum
        W,                  // neutral tungsten
        Re,                 // neutral rhenium
        Os,                 // neutral osmium
        Ir,                 // neutral iridium
        Pt,                 // neutral platinum
        Au,                 // neutral gold
        Tl,                 // neutral thallium
        Pb,                 // neutral lead
        Bi,                 // neutral bismuth
        OTHER,              // all other atoms
        COUNT,              // this will have the numerical value of the number of form factor types, and can thus be used to allocate arrays
        UNKNOWN,            // this is used to indicate that the form factor is unknown
    };

    /**
     * @brief The number of defined form factor types (including excluded volume).
     */
    constexpr int total_ff_count = static_cast<int>(form_factor::form_factor_t::COUNT);

    constexpr int exv_bin   = static_cast<int>(form_factor::form_factor_t::EXCLUDED_VOLUME);
    constexpr int water_bin = static_cast<int>(form_factor::form_factor_t::WATER);
    static_assert(exv_bin == 0, "form_factor::form_factor_t::EXCLUDED_VOLUME must be at index 0");
    static_assert(water_bin == 1, "form_factor::form_factor_t::WATER must be at index 1");

    constexpr int start_index_for_explicit_exv() {return static_cast<int>(form_factor::form_factor_t::EXCLUDED_VOLUME)+1;}

    namespace detail {
        /**
         * @brief Check if a form factor type has a row in the descriptor table, i.e. it is neither COUNT nor UNKNOWN.
         */
        constexpr bool is_tabulated(form_factor_t type) {
            return 0 <= static_cast<int>(type) && static_cast<int>(type) < total_ff_count;
        }
    }
}
