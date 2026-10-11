// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

namespace ausaxs::form_factor {
    // The form factor type of an atom. This is intended to be used as an index for best performance.
    // NOLINTNEXTLINE(readability-enum-initial-value)
    enum class form_factor_t : int {
        EXCLUDED_VOLUME,    // excluded volume
        WATER,              // water
        OH = WATER,         // neutral oxygen with hydrogen
        H,                  // neutral hydrogen
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
        F,                  // neutral fluorine
        Na,                 // neutral sodium
        Mg,                 // neutral magnesium
        P,                  // neutral phosphorus
        Cl,                 // neutral chlorine
        K,                  // neutral potassium
        Ca,                 // neutral calcium
        Mn,                 // neutral manganese
        Fe,                 // neutral iron
        Co,                 // neutral cobalt
        Ni,                 // neutral nickel
        Cu,                 // neutral copper
        Zn,                 // neutral zinc
        Se,                 // neutral selenium
        Br,                 // neutral bromine
        I,                  // neutral iodine
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
