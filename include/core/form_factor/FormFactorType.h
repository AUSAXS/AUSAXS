// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <constants/ConstantsFwd.h>

#include <string>
#include <stdexcept>

namespace ausaxs::form_factor {
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
        NH3,                // neutral nitrogen with three hydrogens
        O,                  // neutral oxygen
        OH,                 // neutral oxygen with hydrogen
        S,                  // neutral sulfur
        SH,                 // neutral sulfur with hydrogen
        OTHER,              // all other atoms
        EXCLUDED_VOLUME,    // excluded volume
        COUNT,              // this will have the numerical value of the number of form factor types, and can thus be used to allocate arrays
        UNKNOWN,            // this is used to indicate that the form factor is unknown
    };
    constexpr int exv_bin   = static_cast<int>(form_factor::form_factor_t::EXCLUDED_VOLUME);
    constexpr int water_bin = static_cast<int>(form_factor::form_factor_t::OH);

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
            case form_factor_t::NH3: return "NH3";
            case form_factor_t::O: return "O";
            case form_factor_t::OH: return "OH";
            case form_factor_t::S: return "S";
            case form_factor_t::SH: return "SH";
            case form_factor_t::OTHER: return "OTH";
            case form_factor_t::EXCLUDED_VOLUME: return "EXV";
            case form_factor_t::COUNT: return "CNT";
            case form_factor_t::UNKNOWN: return "UNK";
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
            case constants::atom_t::S: return form_factor_t::S;
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
            case constants::atomic_group_t::NH3: return form_factor_t::NH3;
            case constants::atomic_group_t::OH: return form_factor_t::OH;
            case constants::atomic_group_t::SH: return form_factor_t::SH;
            default: return get_type(atom_type);
        }
    }

    constexpr constants::atom_t to_atom_type(form_factor_t ff_type) {
        switch(ff_type) {
            case form_factor_t::H: return constants::atom_t::H;
            case form_factor_t::C:
            case form_factor_t::CH:
            case form_factor_t::CH2:
            case form_factor_t::CH3: return constants::atom_t::C;
            case form_factor_t::N:
            case form_factor_t::NH:
            case form_factor_t::NH2:
            case form_factor_t::NH3: return constants::atom_t::N;
            case form_factor_t::O:
            case form_factor_t::OH: return constants::atom_t::O;
            case form_factor_t::S:
            case form_factor_t::SH: return constants::atom_t::S;
            default: return constants::atom_t::Ar;
        }
    }
}

#include <constants/Constants.h>
namespace ausaxs::constants::mass {
    /**
    * @brief Get the mass of an atom in amu.
    */
    constexpr double get_mass(form_factor::form_factor_t type) {
        switch(type) {
            case form_factor::form_factor_t::H: return get_mass(atom_t::H);
            case form_factor::form_factor_t::C: return get_mass(atom_t::C);
            case form_factor::form_factor_t::CH: return 13.019;
            case form_factor::form_factor_t::CH2: return 14.027;
            case form_factor::form_factor_t::CH3: return 15.035;
            case form_factor::form_factor_t::N: return 14.00674;
            case form_factor::form_factor_t::NH: return 15.01474;
            case form_factor::form_factor_t::NH2: return 16.02274;
            case form_factor::form_factor_t::NH3: return 17.03074;
            case form_factor::form_factor_t::O: return 15.999;
            case form_factor::form_factor_t::OH: return 16.999;
            case form_factor::form_factor_t::S: return 32.06;
            case form_factor::form_factor_t::SH: return 33.06;
            case form_factor::form_factor_t::OTHER: return 39.948;
            case form_factor::form_factor_t::EXCLUDED_VOLUME: return 0;
            case form_factor::form_factor_t::COUNT: return 0;
            default: throw std::runtime_error("constants::mass::get_mass: Unknown form factor type \"" + form_factor::to_string(type) + "\"");
        }
    }
}

namespace ausaxs::constants::radius {
    constexpr double get_vdw_radius(form_factor::form_factor_t type) {
        switch(type) {
            case form_factor::form_factor_t::H: return get_vdw_radius(atom_t::H);
            case form_factor::form_factor_t::C: return get_vdw_radius(atom_t::C);
            case form_factor::form_factor_t::CH: return get_vdw_radius(atom_t::C);
            case form_factor::form_factor_t::CH2: return get_vdw_radius(atom_t::C);
            case form_factor::form_factor_t::CH3: return get_vdw_radius(atom_t::C);
            case form_factor::form_factor_t::N: return get_vdw_radius(atom_t::N);
            case form_factor::form_factor_t::NH: return get_vdw_radius(atom_t::N);
            case form_factor::form_factor_t::NH2: return get_vdw_radius(atom_t::N);
            case form_factor::form_factor_t::NH3: return get_vdw_radius(atom_t::N);
            case form_factor::form_factor_t::O: return get_vdw_radius(atom_t::O);
            case form_factor::form_factor_t::OH: return get_vdw_radius(atom_t::O);
            case form_factor::form_factor_t::S: return get_vdw_radius(atom_t::S);
            case form_factor::form_factor_t::SH: return get_vdw_radius(atom_t::S);
            case form_factor::form_factor_t::OTHER: return get_vdw_radius(atom_t::Ar);
            case form_factor::form_factor_t::UNKNOWN: return 0;
            default: throw std::runtime_error("constants::radius::get_vdw_radius: Unknown form factor type \"" + form_factor::to_string(type) + "\"");
        }
    }
}

namespace ausaxs::constants::charge::nuclear {
    /**
     * @brief Get the charge of an atom in e.
     */
    constexpr unsigned int get_charge(form_factor::form_factor_t type) {
        switch(type) {
            case form_factor::form_factor_t::H: return 1;
            case form_factor::form_factor_t::C: return 6;
            case form_factor::form_factor_t::CH: return 7;
            case form_factor::form_factor_t::CH2: return 8;
            case form_factor::form_factor_t::CH3: return 9;
            case form_factor::form_factor_t::N: return 7;
            case form_factor::form_factor_t::NH: return 8;
            case form_factor::form_factor_t::NH2: return 9;
            case form_factor::form_factor_t::NH3: return 10;
            case form_factor::form_factor_t::O: return 8;
            case form_factor::form_factor_t::OH: return 9;
            case form_factor::form_factor_t::S: return 16;
            case form_factor::form_factor_t::SH: return 17;
            case form_factor::form_factor_t::OTHER: return 18;
            case ausaxs::form_factor::form_factor_t::EXCLUDED_VOLUME: return 0;
            default: throw std::runtime_error("constants::charge::nuclear::get_charge: Unknown form factor type \"" + form_factor::to_string(type) + "\"");
        }
    }
}