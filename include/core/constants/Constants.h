// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <constants/ConstantsAxes.h>
#include <constants/ConstantsCoordinates.h>
#include <constants/ConstantsFitParameters.h>
#include <constants/ConstantsFwd.h>
#include <constants/ConstantsProperties.h>
#include <constants/ConstantsSI.h>
#include <constants/ValidFileExtensions.h>
#include <constants/Version.h>
#include <io/IOFwd.h>
#include <math/ConstexprMath.h>
#include <residue/ResidueStorage.h>
#include <utility/Exceptions.h>
#include <utility/SimpleMap.h>

#include <string>

/**
 * @brief This namespace contains all constants used in this project. 
 */
namespace ausaxs::constants {
    /**
     * @brief Radius
     * 
     * This namespace contains all the radius constants used in this project. 
     */
    namespace radius {
        constexpr double electron = 0.0000281794; // electron radius in units of Ångström
    }
    constexpr double Avogadro = 6.02214076e23; // mol^-1

    /**
     * @brief Relative units.
     * 
     * This namespace contains all the unit conversion constants used in this project. 
     * The basic units are for:
     *     mass: Dalton
     *     length: Å
     *     charge: e
     */
    namespace unit {
        constexpr double gm = 1.66054e-24; // Dalton --> grams
        constexpr double mg = 1.66054e-21; // Dalton --> mg
        constexpr double cm = 1e-8; // Ångström --> cm
        constexpr double nm = 1e-1; // Ångström --> nm

        constexpr double mL = math::pow(unit::cm, 3); // Ångström^3 --> mL
    }

    // The 1-symbol names of all amino acids. 
    extern const saxs::detail::SimpleMap<char> name_1symbol_map;

    // The 3-symbol names of all amino acids. 
    extern const saxs::detail::SimpleMap<std::string> name_3symbol_map;

    /**
     * @brief Volume
     * 
     * This namespace contains the volume of all amino acids. 
     * They are taken from https://doi.org/10.1088/0034-4885/39/10/001.
     * All values are in Å^3
     */
    namespace volume {
        // get the volume of a 3symbol amino acid
        extern const saxs::detail::SimpleMap<double> amino_acids;
    }

    /**
     * @brief Radius
     * 
     * This namespace contains the radius of the most common atomic elements encountered in SAXS. 
     */
    namespace radius {
        /**
         * @brief Get the Van der Waals radius of an atom in Ångström.
         */
        double get_vdw_radius(atom_t atom);

        /**
         * @brief Set the radius of the dummy atom.
         */
        void set_dummy_radius(double radius);
        namespace detail {
            extern double dummy_radius;
        }

        constexpr double average_atomic_radius = 1.62;
    }

    /**
     * @brief Charge
     * 
     * This namespace contains the net charge of the most common atomic elements encountered in SAXS. 
     */
    namespace charge {
        namespace nuclear {
            /**
            * @brief Get the charge of an atom in e.
            */
            constexpr int get_charge(atom_t atom);
        }

        namespace ionic {
            /**
             * @brief Get the ionic charge of an atom in e.
             */
            constexpr int get_charge(atom_t atom);
        }

        namespace density {
            constexpr double water = 0.334; // e/Å^3
        }
    }

    /**
     * @brief Valence atoms
     */
    namespace valence {
        /**
         * @brief Get the valence of an atom.
         *        Note that this is *not* the number of valence electrons, but the typical number of bonds the atom can form.
         *        This information is primarily used to determine the number of hydrogens attached to an atom, as these are often not present in the PDB file.
         */
        int get_valence(atom_t atom);
    }

    namespace symbols {
        extern std::string hydrogen;
        extern std::string carbon;
        extern std::string nitrogen;
        extern std::string oxygen;

        namespace detail {
            extern const saxs::detail::SimpleMap<constants::atom_t> string_to_atomt_map;
        }

        atom_t parse_element_string(const std::string& element_string);
        std::string to_string(atom_t atom);
        std::string to_string(atomic_group_t group);

        /**
         * @brief Get the atomic group of an atom.
         * 
         * @param residue_name The name of the residue, e.g. GLY or ALA.
         * @param atom_name The name of the atom, e.g. CH2 or NH2
         * @param atom_type The atom type. This is required to avoid ambiguities since the name is always capitalized in PDB files, so otherwise we cannot distinguish between e.g. a C-alpha (CA) and a calcium (Ca).
         */
        atomic_group_t get_atomic_group(const std::string& residue_name, const std::string& atom_name, atom_t atom_type);

        /**
         * @brief Get the atomic group a given atom belongs to.
         * 
         * @param atom_type The atomic type. 
         * @param hydrogens The number of hydrogens attached to the atom.
         */
        constants::atomic_group_t get_atomic_group(constants::atom_t atom_type, int hydrogens);
    }

    /**
     * @brief Hydrogen atoms
     */
    namespace hydrogen_atoms {
        // get the number of hydrogen atoms attached to an atom of a specific acid. Example: get.at("GLY").at("CA") = 2
        extern residue::ResidueStorage residues;
    }
}

constexpr int ausaxs::constants::charge::nuclear::get_charge(atom_t atom) {
    switch(atom) {
        case atom_t::H:  return 1;
        case atom_t::D:  return 1;
        case atom_t::He: return 2;

        case atom_t::Li: return 3;
        case atom_t::Be: return 4;
        case atom_t::B:  return 5;
        case atom_t::C:  return 6;
        case atom_t::N:  return 7;
        case atom_t::O:  return 8;
        case atom_t::F:  return 9;
        case atom_t::Ne: return 10;

        case atom_t::Na: return 11;
        case atom_t::Mg: return 12;
        case atom_t::Al: return 13;
        case atom_t::Si: return 14;
        case atom_t::P:  return 15;
        case atom_t::S:  return 16;
        case atom_t::Cl: return 17;
        case atom_t::Ar: return 18;

        case atom_t::K:  return 19;
        case atom_t::Ca: return 20;
        case atom_t::Sc: return 21;
        case atom_t::Ti: return 22;
        case atom_t::V:  return 23;
        case atom_t::Cr: return 24;
        case atom_t::Mn: return 25;
        case atom_t::Fe: return 26;
        case atom_t::Co: return 27;
        case atom_t::Ni: return 28;
        case atom_t::Cu: return 29;
        case atom_t::Zn: return 30;
        case atom_t::Ga: return 31;
        case atom_t::Ge: return 32;
        case atom_t::As: return 33;
        case atom_t::Se: return 34;
        case atom_t::Br: return 35;
        case atom_t::Kr: return 36;

        case atom_t::Rb: return 37;
        case atom_t::Sr: return 38;
        case atom_t::Y:  return 39;
        case atom_t::Zr: return 40;
        case atom_t::Nb: return 41;
        case atom_t::Mo: return 42;
        case atom_t::Tc: return 43;
        case atom_t::Ru: return 44;
        case atom_t::Rh: return 45;
        case atom_t::Pd: return 46;
        case atom_t::Ag: return 47;
        case atom_t::Cd: return 48;
        case atom_t::In: return 49;
        case atom_t::Sn: return 50;
        case atom_t::Sb: return 51;
        case atom_t::Te: return 52;
        case atom_t::I:  return 53;
        case atom_t::Xe: return 54;

        case atom_t::Cs: return 55;
        case atom_t::Ba: return 56;
        case atom_t::Hf: return 72;
        case atom_t::Ta: return 73;
        case atom_t::W:  return 74;
        case atom_t::Re: return 75;
        case atom_t::Os: return 76;
        case atom_t::Ir: return 77;
        case atom_t::Pt: return 78;
        case atom_t::Au: return 79;
        case atom_t::Tl: return 81;
        case atom_t::Pb: return 82;
        case atom_t::Bi: return 83;
        case atom_t::Po: return 84;

        case atom_t::M:  return 0;
        case atom_t::dummy: return 1;
        case atom_t::unknown: return 1;
        default: throw ausaxs::except::runtime_error("constants::charge::nuclear::get_charge: Unknown atom type \"" + ausaxs::constants::symbols::to_string(atom) + "\"");
    }
}

constexpr int ausaxs::constants::charge::ionic::get_charge(atom_t atom) {
    switch (atom) {
        case atom_t::Ca: return 2;
        case atom_t::Cl: return -1;
        case atom_t::Zn: return 2;
        default: throw ausaxs::except::runtime_error("constants::charge::ionic::get_charge: Unknown atom type \"" + ausaxs::constants::symbols::to_string(atom) + "\"");
    }
}
