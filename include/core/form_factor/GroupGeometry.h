// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <math/ConstexprMath.h>

#include <numbers>

// The geometry of the united atomic groups, shared by all probes. See each entry for its source.
namespace ausaxs::form_factor::group {
    /**
     * @brief A heavy atom X with n equivalent bound hydrogens: every hydrogen is at the same distance d_xh from X, and every pair of hydrogens is d_hh apart.
     */
    struct Geometry {
        int hydrogens = 0;      // The number of bound hydrogens.
        double d_xh = 0;        // The X-H distance in Å.
        double d_hh = 0;        // The H-H distance in Å. Only meaningful for two or more hydrogens.
    };

    /**
     * @brief X-H distances in Å between the nuclear positions.
     *        These are the same for every probe; a displacement of the hydrogen electron into the bond belongs to the X-ray form factor of the hydrogen.
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

    constexpr Geometry atom {}; // a lone atom without hydrogens
    constexpr Geometry CH  {.hydrogens = 1, .d_xh = bond::CH};
    constexpr Geometry CH2 {.hydrogens = 2, .d_xh = bond::CH, .d_hh = hh_distance(bond::CH, angle::tetrahedral)};
    constexpr Geometry CH3 {.hydrogens = 3, .d_xh = bond::CH, .d_hh = hh_distance(bond::CH, angle::tetrahedral)};
    constexpr Geometry NH  {.hydrogens = 1, .d_xh = bond::NH};
    constexpr Geometry NH2 {.hydrogens = 2, .d_xh = bond::NH, .d_hh = hh_distance(bond::NH, angle::trigonal)}; // only the planar amide and guanidinium groups of proteins
    constexpr Geometry NH3 {.hydrogens = 3, .d_xh = bond::NH, .d_hh = hh_distance(bond::NH, angle::tetrahedral)};
    constexpr Geometry OH  {.hydrogens = 1, .d_xh = bond::OH};
    constexpr Geometry SH  {.hydrogens = 1, .d_xh = bond::SH};
}
