#pragma once

#include <constants/SI.h>
#include <constants/vdwTable.h>
#include <math/ConstexprMath.h>
#include <constants/ConstantsMath.h>

namespace ausaxs::constants::displaced_volume {
    namespace detail {
        struct DisplacedVolumeSet {
            double H;
            double C, CH, CH2, CH3;
            double N, NH, NH2, NH3;
            double O, OH;
            double S, SH;
        };

        constexpr double volume(double radius) {
            return 4*constants::pi/3*constexpr_math::pow(radius, 3);
        }
    }

    // from original CRYSOL paper, 1995: https://doi.org/10.1107/S0021889895007047
    constexpr detail::DisplacedVolumeSet Traube {
        .H   = 0.00515*constexpr_math::pow(constants::SI::length::nm/constants::SI::length::A, 3),
        .C   = 0.01644*constexpr_math::pow(constants::SI::length::nm/constants::SI::length::A, 3), 
        .CH  = 0.02159*constexpr_math::pow(constants::SI::length::nm/constants::SI::length::A, 3), 
        .CH2 = 0.02674*constexpr_math::pow(constants::SI::length::nm/constants::SI::length::A, 3), 
        .CH3 = 0.03189*constexpr_math::pow(constants::SI::length::nm/constants::SI::length::A, 3), 
        .N   = 0.00249*constexpr_math::pow(constants::SI::length::nm/constants::SI::length::A, 3), 
        .NH  = 0.00764*constexpr_math::pow(constants::SI::length::nm/constants::SI::length::A, 3), 
        .NH2 = 0.01279*constexpr_math::pow(constants::SI::length::nm/constants::SI::length::A, 3), 
        .NH3 = 0.01794*constexpr_math::pow(constants::SI::length::nm/constants::SI::length::A, 3), 
        .O   = 0.00913*constexpr_math::pow(constants::SI::length::nm/constants::SI::length::A, 3), 
        .OH  = 0.01428*constexpr_math::pow(constants::SI::length::nm/constants::SI::length::A, 3), 
        .S   = 0.01986*constexpr_math::pow(constants::SI::length::nm/constants::SI::length::A, 3), 
        .SH  = 0.02510*constexpr_math::pow(constants::SI::length::nm/constants::SI::length::A, 3)
    };

    // table I, V^vor from Schaefer et al, 2001: https://doi.org/10.1002/JCC.1137
    constexpr detail::DisplacedVolumeSet Voronoi_implicit_H {
        .H   = 0,
        .C   = 8.895,
        .CH  = 12.430,
        .CH2 = 22.033,
        .CH3 = 34.092,
        .N   = 9.558,
        .NH  = 14.944,
        .NH2 = 22.129,
        .NH3 = 20.641,
        .O   = 22.315,
        .OH  = 23.266,
        .S   = 26.356,
        .SH  = 34.192
    };

    // table I, V^mf from Schaefer et al, 2001: https://doi.org/10.1002/JCC.1137
    constexpr detail::DisplacedVolumeSet MinimumFluctuation_implicit_H {
        .H   = 0,
        .C   = 12.352,
        .CH  = 11.640,
        .CH2 = 34.583,
        .CH3 = 41.851,
        .N   = 0.027,
        .NH  = 2.181,
        .NH2 = 20.562,
        .NH3 = 20.722,
        .O   = 14.238,
        .OH  = 20.911,
        .S   = 15.413,
        .SH  = 28.529
    };

    // table I, V^vor_H from Schaefer et al, 2001: https://doi.org/10.1002/JCC.1137
    constexpr detail::DisplacedVolumeSet Voronoi_explicit_H {
        .H   = 12.958,
        .C   = 8.658,
        .CH  = 11.784,
        .CH2 = 20.682,
        .CH3 = 33.175,
        .N   = 9.144,
        .NH  = 7.119,
        .NH2 = 5.859,
        .NH3 = 2.588,
        .O   = 19.167,
        .OH  = 13.099,
        .S   = 25.715,
        .SH  = 32.333
    };

    // table I, V^mf_H from Schaefer et al, 2001: https://doi.org/10.1002/JCC.1137
    constexpr detail::DisplacedVolumeSet MinimumFluctuation_explicit_H {
        .H   = 0.347,
        .C   = 12.734,
        .CH  = 11.399,
        .CH2 = 34.828,
        .CH3 = 42.011,
        .N   = 0.018,
        .NH  = 1.451,
        .NH2 = 19.064,
        .NH3 = 17.498,
        .O   = 14.334,
        .OH  = 20.312,
        .S   = 15.242,
        .SH  = 28.475
    };

    // based on the van der waals radii of each atom
    constexpr detail::DisplacedVolumeSet vdw {
        .H   = detail::volume(constants::radius::vdw::H),
        .C   = detail::volume(constants::radius::vdw::C),
        .CH  = detail::volume(constants::radius::vdw::C) + 1*detail::volume(constants::radius::vdw::H),
        .CH2 = detail::volume(constants::radius::vdw::C) + 2*detail::volume(constants::radius::vdw::H),
        .CH3 = detail::volume(constants::radius::vdw::C) + 3*detail::volume(constants::radius::vdw::H),
        .N   = detail::volume(constants::radius::vdw::N),
        .NH  = detail::volume(constants::radius::vdw::N) + 1*detail::volume(constants::radius::vdw::H),
        .NH2 = detail::volume(constants::radius::vdw::N) + 2*detail::volume(constants::radius::vdw::H),
        .NH3 = detail::volume(constants::radius::vdw::N) + 3*detail::volume(constants::radius::vdw::H),
        .O   = detail::volume(constants::radius::vdw::O),
        .OH  = detail::volume(constants::radius::vdw::O) + detail::volume(constants::radius::vdw::H),
        .S   = detail::volume(constants::radius::vdw::S),
        .SH  = detail::volume(constants::radius::vdw::S) + detail::volume(constants::radius::vdw::H)
    };

    /**
     * @brief Get the currently used displaced volume set as specified by the settings.
     */
    detail::DisplacedVolumeSet get_displaced_volume_set();

    //! Remember to update settings::molecule::DisplacedVolumeSet::Default if this is changed
    inline constexpr const detail::DisplacedVolumeSet& standard = MinimumFluctuation_implicit_H;
    constexpr double OH2 = 2.98*constexpr_math::pow(10., -23)*constexpr_math::pow(constants::SI::length::cm/constants::SI::length::A, 3);
    constexpr double Ar = detail::volume(constants::radius::vdw::Ar);
}