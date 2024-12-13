#pragma once

#include <utility/SimpleMap.h>
#include <residue/ResidueStorage.h>
#include <constants/ConstantsFwd.h>
#include <constants/ConstantsAxes.h>
#include <constants/ConstantsCoords.h>
#include <constants/ConstantsFitParameters.h>
#include <constants/ValidFileExtensions.h>
#include <constants/Version.h>
#include <constants/SI.h>
#include <io/IOFwd.h>
#include <math/ConstexprMath.h>

#include <string>
#include <stdexcept>

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
     * @brief Mass 
     * 
     * This namespace contains the masses of the most common atomic elements encountered in SAXS. 
     */
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
        // get the valence of an atom
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
        constexpr std::string to_string(atom_t atom);
        constexpr std::string to_string(atomic_group_t group);

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
        case atom_t::unknown: throw std::runtime_error("constants::mass::get_mass: Attempting to get mass of \"unknown\" atom type");
        default: throw std::runtime_error("constants::mass:get_mass: Missing switch case for atom type \"" + constants::symbols::to_string(atom) + "\"");
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
        case atomic_group_t::unknown: throw std::runtime_error("constants::mass::get_mass: Attempting to get mass of \"unknown\" atomic group");
        default: throw std::runtime_error("constants::mass::get_mass: Missing switch case for atomic group \"" + constants::symbols::to_string(group) + "\"");
    }
}

constexpr std::string ausaxs::constants::symbols::to_string(atom_t atom) {
    switch(atom) {
        case atom_t::H: return "H";
        case atom_t::He: return "He";
        case atom_t::Li: return "Li";
        case atom_t::Be: return "Be";
        case atom_t::B: return "B";
        case atom_t::C: return "C";
        case atom_t::N: return "N";
        case atom_t::O: return "O";
        case atom_t::F: return "F";
        case atom_t::Ne: return "Ne";
        case atom_t::Na: return "Na";
        case atom_t::Mg: return "Mg";
        case atom_t::Al: return "Al";
        case atom_t::Si: return "Si";
        case atom_t::P: return "P";
        case atom_t::S: return "S";
        case atom_t::Cl: return "Cl";
        case atom_t::Ar: return "Ar";
        case atom_t::K: return "K";
        case atom_t::Ca: return "Ca";
        case atom_t::Sc: return "Sc";
        case atom_t::Ti: return "Ti";
        case atom_t::V: return "V";
        case atom_t::Cr: return "Cr";
        case atom_t::Mn: return "Mn";
        case atom_t::Fe: return "Fe";
        case atom_t::Co: return "Co";
        case atom_t::Ni: return "Ni";
        case atom_t::Cu: return "Cu";
        case atom_t::Zn: return "Zn";
        case atom_t::I: return "I";
        case atom_t::W: return "W";
        case atom_t::M: return "M";
        case atom_t::dummy: return "#";
        default: throw std::runtime_error("constants::symbols::to_string: Unknown atom type \"" + std::to_string(static_cast<int>(atom)) + "\"");
    }
}

constexpr std::string ausaxs::constants::symbols::to_string(atomic_group_t group) {
    switch(group) {
        case atomic_group_t::CH: return "CH";
        case atomic_group_t::CH2: return "CH2";
        case atomic_group_t::CH3: return "CH3";
        case atomic_group_t::NH: return "NH";
        case atomic_group_t::NH2: return "NH2";
        case atomic_group_t::NH3: return "NH3";
        case atomic_group_t::OH: return "OH";
        case atomic_group_t::SH: return "SH";
        case atomic_group_t::unknown: return "unknown";
        default: throw std::runtime_error("constants::symbols::to_string: Unknown atomic group \"" + std::to_string(static_cast<int>(group)) + "\"");
    }
}