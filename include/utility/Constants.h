#pragma once

#include <utility/ResidueParser.h>
#include <utility/SimpleMap.h>

#include <string>
#include <cmath>

namespace io {class ExistingFile;}

/**
 * @brief Constexpr power function.
 */
constexpr double simple_pow(double val, unsigned int power) {
    double sum = 1;
    for (unsigned int i = 0; i < power; i++) {
        sum *= val;
    }
    return sum;
}

/**
 * @brief This namespace contains all constants used in this project. 
 */
namespace constants {
    namespace filetypes {
        namespace detail {
            struct FileType {
                FileType(std::vector<std::string> extensions);
                bool validate(const io::ExistingFile& path) const;
                std::vector<std::string> extensions;
            };
        }

        const detail::FileType structure =  {{".pdb",  ".ent"}};
        const detail::FileType saxs_data =  {{".dat",  ".txt",  ".rsr",  ".xvg"}};
        const detail::FileType em_map =     {{".map",  ".ccp4", ".mrc", ".rec"}};
        const detail::FileType unit_cell =  {{".cell", ".uc"}};
        const detail::FileType grid =       {{".grid"}};
    }

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

        constexpr double mL = simple_pow(unit::cm, 3); // Ångström^3 --> mL
    }

    /**
     * @brief Absolute units.
     * 
     * This namespace contains all the absolute unit conversion constants. 
     */
    namespace SI {
        namespace mass {
            constexpr double kg = 1;
            constexpr double gm = 1e-3;
            constexpr double mg = 1e-6;
            constexpr double u = 1.66053*1e-27;
        }

        namespace length {
            constexpr double m = 1;
            constexpr double cm = 1e-2;
            constexpr double A = 1e-10; // Ångström
        }

        namespace volume {
            constexpr double A3 = 1e-30; // Ångström^3
            constexpr double cm3 = 1e-6;
        }
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

    // Atom enum for a consistent way to denote the type of an atom
    enum class atom_t {
        e, H, He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, P, S, Cl, Ar, K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn, W, 
        M,      // This fake element is for compatibility with the GROMACS tip4p water model.
        dummy,  // Fake element for dummy atoms with variable properties. 
        unknown
    };

    /**
     * @brief Mass 
     * 
     * This namespace contains the masses of the most common atomic elements encountered in SAXS. 
     */
    namespace mass {
        /**
         * @brief Get the mass of an atom in u.
         */
        double get_mass(atom_t atom);

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
    }

    /**
     * @brief Charge
     * 
     * This namespace contains the net charge of the most common atomic elements encountered in SAXS. 
     */
    namespace charge {
        /**
         * @brief Get the charge of an atom in e.
         */
        unsigned int get_charge(atom_t atom);

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
        std::string write_element_string(atom_t atom);
    }

    /**
     * @brief Hydrogen atoms
     */
    namespace hydrogen_atoms {
        // get the number of hydrogen atoms attached to an atom of a specific acid. Example: get.at("GLY").at("CA") = 2
        extern parser::residue::ResidueStorage residues;
    }

    namespace form_factor {
        constexpr double sigma_excluded_volume = 1.62; // Å
        constexpr double fpi2 = 4*M_PI*M_PI; // 4π^2, conversion factor from s to q

        // no source
        namespace hydrogen {
            constexpr std::array<double, 5> a = {1, 0, 0, 0, 0};
            constexpr std::array<double, 5> b = {0.5, 0, 0, 0, 0};
            constexpr double c = 0;
        }

        // Waasmeier & Kirfel, https://doi.org/10.1107/S0108767394013292
        // this is the form factor of argon
        namespace other {
            constexpr std::array<double, 5> a = {7.188004, 6.638454, 0.454180, 1.929593, 1.523654};
            constexpr std::array<double, 5> b = {0.956221/fpi2, 15.339877/fpi2, 15.339862/fpi2, 39.043824/fpi2, 0.062409/fpi2};
            constexpr double c = 0.265954;
        }

        //! get source from Jan
        namespace excluded_volume {
            constexpr std::array<double, 5> a = {1, 0, 0, 0, 0};
            constexpr std::array<double, 5> b = {1.62*1.62/2, 0, 0, 0, 0};
            constexpr double c = 0;
        }

        // Waasmeier & Kirfel, https://doi.org/10.1107/S0108767394013292
        namespace neutral_carbon {
            constexpr std::array<double, 5> a = { 2.657506, 1.078079,  1.490909, -4.241070, 0.713791};
            constexpr std::array<double, 5> b = {14.780758/fpi2, 0.776775/fpi2, 42.086843/fpi2, -0.000294/fpi2, 0.239535/fpi2};
            constexpr double c = 4.297983;
        }

        // Waasmeier & Kirfel, https://doi.org/10.1107/S0108767394013292
        namespace neutral_oxygen {
            constexpr std::array<double, 5> a = { 2.960427, 2.508818, 0.637853,  0.722838, 1.142756};
            constexpr std::array<double, 5> b = {14.182259/fpi2, 5.936858/fpi2, 0.112726/fpi2, 34.958481/fpi2, 0.390240/fpi2};
            constexpr double c = 0.027014;
        }

        // Waasmeier & Kirfel, https://doi.org/10.1107/S0108767394013292
        namespace neutral_nitrogen {
            constexpr std::array<double, 5> a = {11.893780,  3.277479,  1.858092, 0.858927, 0.912985};
            constexpr std::array<double, 5> b = { 0.000158/fpi2, 10.232723/fpi2, 30.344690/fpi2, 0.656065/fpi2, 0.217287/fpi2};
            constexpr double c = -11.804902;
        }
    }
}