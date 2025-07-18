// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <utility/SimpleMap.h>
#include <residue/ResidueStorage.h>
#include <constants/ConstantsFwd.h>
#include <constants/ConstantsAxes.h>
#include <constants/ConstantsCoordinates.h>
#include <constants/ConstantsFitParameters.h>
#include <constants/ConstantsProperties.h>
#include <constants/ValidFileExtensions.h>
#include <constants/Version.h>
#include <constants/ConstantsSI.h>
#include <io/IOFwd.h>
#include <math/ConstexprMath.h>

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
    constexpr double Avogadro = 6.02214076e-23; // mol^-1

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

        constexpr double mL = constexpr_math::pow(unit::cm, 3); // Ångström^3 --> mL
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
            unsigned int get_charge(atom_t atom);
        }

        namespace ionic {
            /**
             * @brief Get the ionic charge of an atom in e.
             */
            unsigned int get_charge(atom_t atom);
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
        unsigned int get_valence(atom_t atom);
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
        constants::atomic_group_t get_atomic_group(constants::atom_t atom_type, unsigned int hydrogens);
    }

    /**
     * @brief Hydrogen atoms
     */
    namespace hydrogen_atoms {
        // get the number of hydrogen atoms attached to an atom of a specific acid. Example: get.at("GLY").at("CA") = 2
        extern residue::ResidueStorage residues;
    }
}