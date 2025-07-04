// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <constants/ConstantsFwd.h>
#include <constants/ConstantsSI.h>

#include <cassert>
#include <stdexcept>

namespace ausaxs::constants {
    namespace mass {
        /**
         * @brief Get the mass of an atom in u.
         */
        constexpr double get_mass(atom_t atom);

        /**
         * @brief Get the mass of an atomic group in u.
         */
        constexpr double get_mass(atomic_group_t group); 

        namespace density {
            constexpr double water = 0.9982067*SI::mass::u/SI::volume::A3;
            constexpr double protein = 1.35*SI::mass::gm/SI::volume::cm3;
        }
    }
}

constexpr double ausaxs::constants::mass::get_mass(atom_t atom) {
    switch(atom) {
        case atom_t::H: return 1.0079;
        case atom_t::He: return 4.0026;
        case atom_t::Li: return 6.941;
        case atom_t::Be: return 9.0122;
        case atom_t::B: return 10.811;
        case atom_t::C: return 12.0107;
        case atom_t::N: return 14.0067;
        case atom_t::O: return 15.9994;
        case atom_t::F: return 18.9984;
        case atom_t::Ne: return 20.1797;
        case atom_t::Na: return 22.9897;
        case atom_t::Mg: return 24.305;
        case atom_t::Al: return 26.9815;
        case atom_t::Si: return 28.0855;
        case atom_t::P: return 30.9738;
        case atom_t::S: return 32.065;
        case atom_t::Cl: return 35.453;
        case atom_t::Ar: return 39.948;
        case atom_t::K: return 39.0983;
        case atom_t::Ca: return 40.078;
        case atom_t::Sc: return 44.9559;
        case atom_t::Ti: return 47.867;
        case atom_t::V: return 50.9415;
        case atom_t::Cr: return 51.9961;
        case atom_t::Mn: return 54.938;
        case atom_t::Fe: return 55.845;
        case atom_t::Co: return 58.9332;
        case atom_t::Ni: return 58.6934;
        case atom_t::Cu: return 63.546;
        case atom_t::Zn: return 65.39;
        case atom_t::I: return 126.904;
        case atom_t::W: return 183.84;
        case atom_t::M: return 0;
        case atom_t::dummy: return 1;
        case atom_t::unknown: 
            assert(false && "constants::mass::get_mass: Attempting to get mass of \"unknown\" atom type");
            [[fallthrough]];
        default: 
            throw std::runtime_error("constants::mass::get_mass: Missing switch case for atom type: " + std::to_string(static_cast<int>(atom)));
        return 0;
    }
}

constexpr double ausaxs::constants::mass::get_mass(atomic_group_t group) {
    switch(group) {
        case atomic_group_t::CH: return get_mass(atom_t::C) + get_mass(atom_t::H);
        case atomic_group_t::CH2: return get_mass(atom_t::C) + 2*get_mass(atom_t::H);
        case atomic_group_t::CH3: return get_mass(atom_t::C) + 3*get_mass(atom_t::H);
        case atomic_group_t::NH: return get_mass(atom_t::N) + get_mass(atom_t::H);
        case atomic_group_t::NH2: return get_mass(atom_t::N) + 2*get_mass(atom_t::H);
        case atomic_group_t::NH3: return get_mass(atom_t::N) + 3*get_mass(atom_t::H);
        case atomic_group_t::OH: return get_mass(atom_t::O) + get_mass(atom_t::H);
        case atomic_group_t::SH: return get_mass(atom_t::S) + get_mass(atom_t::H);
        case atomic_group_t::unknown: 
            assert(false && "constants::mass::get_mass: Attempting to get mass of \"unknown\" atomic group");
            [[fallthrough]];
        default: 
            throw std::runtime_error("constants::mass::get_mass: Missing switch case for atomic group: " + std::to_string(static_cast<int>(group)));
        return 0;
    }
}