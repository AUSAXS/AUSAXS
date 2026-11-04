// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <math/ConstexprMath.h>

#include <numbers>

// Neutron scattering lengths and group geometries. See each entry for its source.
namespace ausaxs::form_factor::neutron {
    /**
     * @brief Coherent scattering lengths in fm, as distributed by periodictable 2.0.2 (also used by SasView).
     *        Natural elements are from Rauch & Waschkowski, ILL Neutron Data Booklet, 2nd ed. (2003),
     *        with newer measurements for 1H, 2H, C, O, Sn, Pb, Bi from Snow et al., Phys. Rev. D 101, 062004 (2020), https://doi.org/10.1103/PhysRevD.101.062004
     *        and for 4He from Haun et al., Phys. Rev. Lett. 124, 012501 (2020), https://doi.org/10.1103/PhysRevLett.124.012501 (natural He is its isotope-weighted average with 3He).
     *        H is protium, not natural hydrogen, since every hydrogen is treated as protium.
     */
    namespace scattering_length {
        constexpr double H  = -3.7395;
        constexpr double D  =  6.6681;
        constexpr double He =  3.0985;
        constexpr double Li =    -1.9;
        constexpr double Be =    7.79;
        constexpr double B  =     5.3;
        constexpr double C  =  6.6472;
        constexpr double N  =    9.36;
        constexpr double O  =  5.8037;
        constexpr double F  =   5.654;
        constexpr double Ne =   4.566;
        constexpr double Na =    3.63;
        constexpr double Mg =   5.375;
        constexpr double Al =   3.449;
        constexpr double Si = 4.15071;
        constexpr double P  =    5.13;
        constexpr double S  =   2.847;
        constexpr double Cl =  9.5792;
        constexpr double Ar =   1.909;
        constexpr double K  =    3.67;
        constexpr double Ca =     4.7;
        constexpr double Sc =    12.1;
        constexpr double Ti =   -3.37;
        constexpr double V  =  -0.443;
        constexpr double Cr =   3.635;
        constexpr double Mn =   -3.75;
        constexpr double Fe =    9.45;
        constexpr double Co =    2.49;
        constexpr double Ni =    10.3;
        constexpr double Cu =   7.718;
        constexpr double Zn =    5.68;
        constexpr double Ga =   7.288;
        constexpr double Ge =   8.185;
        constexpr double As =    6.58;
        constexpr double Se =    7.97;
        constexpr double Br =    6.79;
        constexpr double Kr =    7.81;
        constexpr double Rb =    7.08;
        constexpr double Sr =    7.02;
        constexpr double Y  =    7.75;
        constexpr double Zr =    7.16;
        constexpr double Nb =   7.054;
        constexpr double Mo =   6.715;
        constexpr double Tc =     6.8;
        constexpr double Ru =    7.02;
        constexpr double Rh =     5.9;
        constexpr double Pd =    5.91;
        constexpr double Ag =   5.922;
        constexpr double Cd =    4.83;
        constexpr double In =   4.065;
        constexpr double Sn =  6.2239;
        constexpr double Sb =    5.57;
        constexpr double Te =    5.68;
        constexpr double I  =    5.28;
        constexpr double Xe =    4.69;
        constexpr double Cs =    5.42;
        constexpr double Ba =    5.07;
        constexpr double Hf =    7.77;
        constexpr double Ta =    6.91;
        constexpr double W  =   4.755;
        constexpr double Re =     9.2;
        constexpr double Os =    10.7;
        constexpr double Ir =    10.6;
        constexpr double Pt =     9.6;
        constexpr double Au =     7.9;
        constexpr double Tl =   8.776;
        constexpr double Pb =  9.4024;
        constexpr double Bi =  8.5242;
    }

    namespace coefficients {
        /**
         * @brief A heavy atom X with n equivalent bound hydrogens: every hydrogen is at the same distance d_xh from X, and every pair of hydrogens is d_hh apart.
         */
        struct Group {
            double b_heavy = 0;     // The coherent scattering length of the heavy atom in fm.
            int hydrogens = 0;      // The number of bound hydrogens.
            double d_xh = 0;        // The X-H distance in Å.
            double d_hh = 0;        // The H-H distance in Å. Only meaningful for two or more hydrogens.
        };

        /**
         * @brief X-H distances in Å between the nuclear positions, which are longer than the X-ray distances since the hydrogen electron is displaced into the bond.
         *        C-H, N-H, and O-H are the neutron-normalized values of Allen, Acta Cryst. B42, 515-522 (1986).
         *        S-H is the neutron value of Allen & Bruno, Acta Cryst. B66, 380-386 (2010).
         */
        namespace bond {
            constexpr double CH = 1.083;
            constexpr double NH = 1.009;
            constexpr double OH = 0.983;
            constexpr double SH = 1.338;
        }

        /**
         * @brief H-X-H angles in degrees.
         */
        namespace angle {
            constexpr double tetrahedral = 109.4712;
            constexpr double trigonal = 120;
        }

        /**
         * @brief The H-H distance of two hydrogens at distance d_xh from X, separated by the H-X-H angle in degrees.
         */
        constexpr double hh_distance(double d_xh, double hh_angle) {
            return 2*d_xh*math::sin(hh_angle*std::numbers::pi/360);
        }

        // Only the atomic groups are defined here. A bare atom needs nothing but its scattering length, e.g. Group{.b_heavy = scattering_length::Ag}.
        namespace b = scattering_length;
        constexpr Group CH  {.b_heavy = b::C, .hydrogens = 1, .d_xh = bond::CH};
        constexpr Group CH2 {.b_heavy = b::C, .hydrogens = 2, .d_xh = bond::CH, .d_hh = hh_distance(bond::CH, angle::tetrahedral)};
        constexpr Group CH3 {.b_heavy = b::C, .hydrogens = 3, .d_xh = bond::CH, .d_hh = hh_distance(bond::CH, angle::tetrahedral)};
        constexpr Group NH  {.b_heavy = b::N, .hydrogens = 1, .d_xh = bond::NH};
        constexpr Group NH2 {.b_heavy = b::N, .hydrogens = 2, .d_xh = bond::NH, .d_hh = hh_distance(bond::NH, angle::trigonal)}; // only the planar amide and guanidinium groups of proteins
        constexpr Group NH3 {.b_heavy = b::N, .hydrogens = 3, .d_xh = bond::NH, .d_hh = hh_distance(bond::NH, angle::tetrahedral)};
        constexpr Group OH  {.b_heavy = b::O, .hydrogens = 1, .d_xh = bond::OH};
        constexpr Group SH  {.b_heavy = b::S, .hydrogens = 1, .d_xh = bond::SH};
    }
}
