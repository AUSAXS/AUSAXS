#pragma once

#include <constants/ConstantsFwd.h>

#include <string>
#include <stdexcept>

namespace form_factor {
    // The form factor type of an atom. This is intended to be used as an index for best performance.
    enum class form_factor_t {
        H,                  // neutral hydrogen
        C,                  // neutral carbon
        CH,                 // neutral carbon with hydrogen
        CH2,                // neutral carbon with two hydrogens
        CH3,                // neutral carbon with three hydrogens
        N,                  // neutral nitrogen
        NH,                 // neutral nitrogen with hydrogen
        NH2,                // neutral nitrogen with two hydrogens
        O,                  // neutral oxygen
        OH,                 // neutral oxygen with hydrogen
        S,                  // neutral sulfur
        SH,                 // neutral sulfur with hydrogen
        OTHER,              // all other atoms
        EXCLUDED_VOLUME,    // excluded volume
        COUNT,              // this will have the numerical value of the number of form factor types, and can thus be used to allocate arrays
    };

    [[maybe_unused]] static std::string to_string(form_factor_t type) {
        switch (type) {
            case form_factor_t::H: return "H";
            case form_factor_t::C: return "C";
            case form_factor_t::CH: return "CH";
            case form_factor_t::CH2: return "CH2";
            case form_factor_t::CH3: return "CH3";
            case form_factor_t::N: return "N";
            case form_factor_t::NH: return "NH";
            case form_factor_t::NH2: return "NH2";
            case form_factor_t::O: return "O";
            case form_factor_t::OH: return "OH";
            case form_factor_t::S: return "S";
            case form_factor_t::SH: return "SH";
            case form_factor_t::OTHER: return "other";
            case form_factor_t::EXCLUDED_VOLUME: return "excluded volume";
            case form_factor_t::COUNT: return "count";
            default: throw std::runtime_error("form_factor::to_string: Invalid form factor type (enum " + std::to_string(static_cast<int>(type)) + ")");
        }
    }

    /**
     * @brief Get the number of unique form factors.
     * 
     * This can be used to iterate over all form factors.
     */
    constexpr unsigned int get_count() {
        return static_cast<unsigned int>(form_factor_t::COUNT);
    }

    /**
     * @brief Get the number of unique form factors, excluding the excluded volume form factor.
     * 
     * This can be used to iterate over all form factors except the excluded volume form factor.
     */
    constexpr unsigned int get_count_without_excluded_volume() {
        return static_cast<unsigned int>(form_factor_t::COUNT)-1;
    }

    /**
     * @brief Get the form factor type based on an atom type. 
     *        In case the atom type is not recognized, the default form factor (argon) is returned.
     */
    constexpr form_factor_t get_type(constants::atom_t atom_type) {
        switch(atom_type) {
            case constants::atom_t::H: return form_factor_t::H;
            case constants::atom_t::C: return form_factor_t::C;
            case constants::atom_t::N: return form_factor_t::N;
            case constants::atom_t::O: return form_factor_t::O;
            default: return form_factor_t::OTHER;
        }
    }

    /**
     * @brief Get the form factor type based on an atom type and an atomic group.
     *        The atomic group takes priority. Only if the atomic group is not recognized, the atom type is used.
     *        In case either the atomic group or the atom type is not recognized, the default form factor (argon) is returned.
     */
    constexpr form_factor_t get_type(constants::atom_t atom_type, constants::atomic_group_t atomic_group) {
        switch(atomic_group) {
            case constants::atomic_group_t::CH: return form_factor_t::CH;
            case constants::atomic_group_t::CH2: return form_factor_t::CH2;
            case constants::atomic_group_t::CH3: return form_factor_t::CH3;
            case constants::atomic_group_t::NH: return form_factor_t::NH;
            case constants::atomic_group_t::NH2: return form_factor_t::NH2;
            case constants::atomic_group_t::OH: return form_factor_t::OH;
            case constants::atomic_group_t::SH: return form_factor_t::SH;
            default: return get_type(atom_type);
        }
    }
}