// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <math/ConstexprMath.h>

#include <numbers>

// Neutron scattering lengths and group geometries. See each entry for its source.
namespace ausaxs::form_factor::neutron {
    /**
     * @brief Coherent scattering lengths in fm.
     *        Sears, Neutron News 3(3), 26-37 (1992), https://doi.org/10.1080/10448639208218770
     */
    namespace scattering_length {
        constexpr double H  = -3.7390;
        constexpr double C  =  6.6460;
        constexpr double N  =  9.36;
        constexpr double O  =  5.803;
        constexpr double F  =  5.654;
        constexpr double Na =  3.63;
        constexpr double Mg =  5.375;
        constexpr double P  =  5.13;
        constexpr double S  =  2.847;
        constexpr double Cl =  9.5770;
        constexpr double Ar =  1.909;
        constexpr double K  =  3.67;
        constexpr double Ca =  4.70;
        constexpr double Mn = -3.73;
        constexpr double Fe =  9.45;
        constexpr double Co =  2.49;
        constexpr double Ni = 10.3;
        constexpr double Cu =  7.718;
        constexpr double Zn =  5.680;
        constexpr double Se =  7.970;
        constexpr double Br =  6.795;
        constexpr double I  =  5.28;
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

        namespace b = scattering_length;
        constexpr Group H   {.b_heavy = b::H};
        constexpr Group C   {.b_heavy = b::C};
        constexpr Group CH  {.b_heavy = b::C, .hydrogens = 1, .d_xh = bond::CH};
        constexpr Group CH2 {.b_heavy = b::C, .hydrogens = 2, .d_xh = bond::CH, .d_hh = hh_distance(bond::CH, angle::tetrahedral)};
        constexpr Group CH3 {.b_heavy = b::C, .hydrogens = 3, .d_xh = bond::CH, .d_hh = hh_distance(bond::CH, angle::tetrahedral)};
        constexpr Group N   {.b_heavy = b::N};
        constexpr Group NH  {.b_heavy = b::N, .hydrogens = 1, .d_xh = bond::NH};
        constexpr Group NH2 {.b_heavy = b::N, .hydrogens = 2, .d_xh = bond::NH, .d_hh = hh_distance(bond::NH, angle::trigonal)}; // only the planar amide and guanidinium groups of proteins
        constexpr Group NH3 {.b_heavy = b::N, .hydrogens = 3, .d_xh = bond::NH, .d_hh = hh_distance(bond::NH, angle::tetrahedral)};
        constexpr Group O   {.b_heavy = b::O};
        constexpr Group OH  {.b_heavy = b::O, .hydrogens = 1, .d_xh = bond::OH};
        constexpr Group S   {.b_heavy = b::S};
        constexpr Group SH  {.b_heavy = b::S, .hydrogens = 1, .d_xh = bond::SH};
        constexpr Group F   {.b_heavy = b::F};
        constexpr Group Na  {.b_heavy = b::Na};
        constexpr Group Mg  {.b_heavy = b::Mg};
        constexpr Group P   {.b_heavy = b::P};
        constexpr Group Cl  {.b_heavy = b::Cl};
        constexpr Group Ar  {.b_heavy = b::Ar};
        constexpr Group K   {.b_heavy = b::K};
        constexpr Group Ca  {.b_heavy = b::Ca};
        constexpr Group Mn  {.b_heavy = b::Mn};
        constexpr Group Fe  {.b_heavy = b::Fe};
        constexpr Group Co  {.b_heavy = b::Co};
        constexpr Group Ni  {.b_heavy = b::Ni};
        constexpr Group Cu  {.b_heavy = b::Cu};
        constexpr Group Zn  {.b_heavy = b::Zn};
        constexpr Group Se  {.b_heavy = b::Se};
        constexpr Group Br  {.b_heavy = b::Br};
        constexpr Group I   {.b_heavy = b::I};
        constexpr Group other = Ar; // default form factor for unknown atoms

        // the excluded volume has no scattering length
        constexpr Group excluded_volume {};
    }
}
