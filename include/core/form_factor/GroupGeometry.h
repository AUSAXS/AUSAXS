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
     * @brief X-H distances in Å between the nuclear positions, one for each group.
     *        These are the same for every probe. The X-ray hydrogen is a spherical density centred on the nucleus (see xray/XrayFormFactorTable.h),
     *        so the displacement of its electron into the bond is not modelled.
     *
     *        The values are the mean neutron-diffraction distances of Allen & Bruno, Acta Cryst. B66, 380-386 (2010), https://doi.org/10.1107/S0108768110012048,
     *        using the most specific bond class matching the group, and the precise subset of their table 4 when it contains that class, otherwise their table 2.
     *        They are transcribed from the supplementary data of Woińska et al., Sci. Adv. 2, e1600192 (2016), https://doi.org/10.6084/m9.figshare.3205588.v1,
     *        which lists the mean, sample standard deviation, and number of observations of every class.
     *        S-H is not covered by those tables and is the neutron value of Allen & Bruno (2010).
     */
    namespace bond {
        constexpr double CH  = 1.100; // Z3-Csp3-H, table 4. Aromatic C(ar)-H is also represented by CH, but is 0.015 Å shorter (1.085, table 4).
        constexpr double CH2 = 1.097; // Z2-Csp3-H2, table 4
        constexpr double CH3 = 1.087; // Z-Csp3-H3, table 4
        constexpr double NH  = 1.027; // Z2-N-H, table 2
        constexpr double NH2 = 1.013; // Csp2-N-H2, table 2, covering both amide and guanidinium groups
        constexpr double NH3 = 1.040; // N+-H, table 4. Our NH3 is always a cation, see form_factor/FormFactorInfo.h.
        constexpr double OH  = 0.989; // Z-O-H, table 4
        constexpr double SH  = 1.338;
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
    constexpr Geometry CH2 {.hydrogens = 2, .d_xh = bond::CH2, .d_hh = hh_distance(bond::CH2, angle::tetrahedral)};
    constexpr Geometry CH3 {.hydrogens = 3, .d_xh = bond::CH3, .d_hh = hh_distance(bond::CH3, angle::tetrahedral)};
    constexpr Geometry NH  {.hydrogens = 1, .d_xh = bond::NH};
    constexpr Geometry NH2 {.hydrogens = 2, .d_xh = bond::NH2, .d_hh = hh_distance(bond::NH2, angle::trigonal)}; // only the planar amide and guanidinium groups of proteins
    constexpr Geometry NH3 {.hydrogens = 3, .d_xh = bond::NH3, .d_hh = hh_distance(bond::NH3, angle::tetrahedral)};
    constexpr Geometry OH  {.hydrogens = 1, .d_xh = bond::OH};
    constexpr Geometry SH  {.hydrogens = 1, .d_xh = bond::SH};
}
