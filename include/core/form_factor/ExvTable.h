// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <constants/ConstantsSI.h>
#include <constants/vdwTable.h>
#include <form_factor/FormFactorType.h>
#include <math/ConstexprMath.h>

#include <array>
#include <numbers>
#include <optional>

// Per-type displaced solvent volumes for the various excluded-volume sets.
// Each set is annotated with its literature source below; see settings::exv::ExvSet for selection.
namespace ausaxs::constants::exv {
    namespace detail {
        /**
         * @brief A set of displaced solvent volumes, with one optional entry per form factor type (a bare atom,
         *        or an atom with its implicit hydrogens, e.g. CH3). Volumes are stored in Å³.
         *        Types without an entry cannot be used with the Fraser excluded volume model.
         */
        struct ExvSet {
            std::array<std::optional<double>, ausaxs::form_factor::total_ff_count> volumes;
            constexpr bool operator==(const ExvSet& other) const = default;

            /**
             * @brief Check if this set has a displaced solvent volume for the given form factor type.
             */
            constexpr bool contains(ausaxs::form_factor::form_factor_t type) const {
                return ausaxs::form_factor::detail::is_tabulated(type) && volumes[static_cast<int>(type)].has_value();
            }

            /**
             * @brief Get the displaced solvent volume of a single atom of the given form factor type, in cubic angstroms.
             */
            constexpr double get(ausaxs::form_factor::form_factor_t type) const {
                if (!contains(type)) {
                    throw ausaxs::except::runtime_error(
                        "constants::exv::detail::ExvSet::get: No displaced volume for form factor type \"" + ausaxs::form_factor::to_string(type) + "\"");
                }
                return *volumes[static_cast<int>(type)];
            }
        };

        /**
         * @brief Descriptor of the displaced solvent volumes of a single form factor type across all volume sets.
         *        Each volume is optional; a missing volume means the type is absent from that set.
         */
        struct ExvInfo {
            ausaxs::form_factor::form_factor_t type;
            std::optional<double> Traube;
            std::optional<double> Voronoi_implicit_H;
            std::optional<double> MinimumFluctuation_implicit_H;
            std::optional<double> Voronoi_explicit_H;
            std::optional<double> MinimumFluctuation_explicit_H;
            std::optional<double> vdw;
        };

        constexpr double volume(double radius) {
            return 4*std::numbers::pi/3*constexpr_math::pow(radius, 3);
        }

        constexpr double nm3(double V) {
            return V*constexpr_math::pow(constants::SI::length::nm/constants::SI::length::A, 3);
        }

        namespace vdw = constants::radius::vdw;
        using ff_t = ausaxs::form_factor::form_factor_t;

        /**
         * @brief The excluded volume descriptor table.
         *        Rows may appear in any order, and form factor types without any known volumes can be omitted entirely.
         *        OTHER must always be present, since it is the fallback type for everything else.
         *
         * Sources of each column:
         *   Traube:                        original CRYSOL paper, 1995: https://doi.org/10.1107/S0021889895007047
         *   Voronoi_implicit_H:            table I, V^vor   from Schaefer et al, 2001: https://doi.org/10.1002/JCC.1137
         *   MinimumFluctuation_implicit_H: table I, V^mf    from Schaefer et al, 2001: https://doi.org/10.1002/JCC.1137
         *   Voronoi_explicit_H:            table I, V^vor_H from Schaefer et al, 2001: https://doi.org/10.1002/JCC.1137
         *   MinimumFluctuation_explicit_H: table I, V^mf_H  from Schaefer et al, 2001: https://doi.org/10.1002/JCC.1137
         *   vdw:                           based on the van der Waals radii of each atom
         */
        constexpr std::array exv_info_table = {
            ExvInfo{.type=ff_t::H,   .Traube=nm3(0.00515), .Voronoi_implicit_H=0,      .MinimumFluctuation_implicit_H=0,      .Voronoi_explicit_H=12.958, .MinimumFluctuation_explicit_H=0.347,  .vdw=volume(vdw::H)                   },
            ExvInfo{.type=ff_t::C,   .Traube=nm3(0.01644), .Voronoi_implicit_H=8.895,  .MinimumFluctuation_implicit_H=12.352, .Voronoi_explicit_H=8.658,  .MinimumFluctuation_explicit_H=12.734, .vdw=volume(vdw::C)                   },
            ExvInfo{.type=ff_t::CH,  .Traube=nm3(0.02159), .Voronoi_implicit_H=12.430, .MinimumFluctuation_implicit_H=11.640, .Voronoi_explicit_H=11.784, .MinimumFluctuation_explicit_H=11.399, .vdw=volume(vdw::C) + 1*volume(vdw::H)},
            ExvInfo{.type=ff_t::CH2, .Traube=nm3(0.02674), .Voronoi_implicit_H=22.033, .MinimumFluctuation_implicit_H=34.583, .Voronoi_explicit_H=20.682, .MinimumFluctuation_explicit_H=34.828, .vdw=volume(vdw::C) + 2*volume(vdw::H)},
            ExvInfo{.type=ff_t::CH3, .Traube=nm3(0.03189), .Voronoi_implicit_H=34.092, .MinimumFluctuation_implicit_H=41.851, .Voronoi_explicit_H=33.175, .MinimumFluctuation_explicit_H=42.011, .vdw=volume(vdw::C) + 3*volume(vdw::H)},
            ExvInfo{.type=ff_t::N,   .Traube=nm3(0.00249), .Voronoi_implicit_H=9.558,  .MinimumFluctuation_implicit_H=0.027,  .Voronoi_explicit_H=9.144,  .MinimumFluctuation_explicit_H=0.018,  .vdw=volume(vdw::N)                   },
            ExvInfo{.type=ff_t::NH,  .Traube=nm3(0.00764), .Voronoi_implicit_H=14.944, .MinimumFluctuation_implicit_H=2.181,  .Voronoi_explicit_H=7.119,  .MinimumFluctuation_explicit_H=1.451,  .vdw=volume(vdw::N) + 1*volume(vdw::H)},
            ExvInfo{.type=ff_t::NH2, .Traube=nm3(0.01279), .Voronoi_implicit_H=22.129, .MinimumFluctuation_implicit_H=20.562, .Voronoi_explicit_H=5.859,  .MinimumFluctuation_explicit_H=19.064, .vdw=volume(vdw::N) + 2*volume(vdw::H)},
            ExvInfo{.type=ff_t::NH3, .Traube=nm3(0.01794), .Voronoi_implicit_H=20.641, .MinimumFluctuation_implicit_H=20.722, .Voronoi_explicit_H=2.588,  .MinimumFluctuation_explicit_H=17.498, .vdw=volume(vdw::N) + 3*volume(vdw::H)},
            ExvInfo{.type=ff_t::O,   .Traube=nm3(0.00913), .Voronoi_implicit_H=22.315, .MinimumFluctuation_implicit_H=14.238, .Voronoi_explicit_H=19.167, .MinimumFluctuation_explicit_H=14.334, .vdw=volume(vdw::O)                   },
            ExvInfo{.type=ff_t::OH,  .Traube=nm3(0.01428), .Voronoi_implicit_H=23.266, .MinimumFluctuation_implicit_H=20.911, .Voronoi_explicit_H=13.099, .MinimumFluctuation_explicit_H=20.312, .vdw=volume(vdw::O) + volume(vdw::H)  },
            ExvInfo{.type=ff_t::S,   .Traube=nm3(0.01986), .Voronoi_implicit_H=26.356, .MinimumFluctuation_implicit_H=15.413, .Voronoi_explicit_H=25.715, .MinimumFluctuation_explicit_H=15.242, .vdw=volume(vdw::S)                   },
            ExvInfo{.type=ff_t::SH,  .Traube=nm3(0.02510), .Voronoi_implicit_H=34.192, .MinimumFluctuation_implicit_H=28.529, .Voronoi_explicit_H=32.333, .MinimumFluctuation_explicit_H=28.475, .vdw=volume(vdw::S) + volume(vdw::H)  },

            // all other atoms are treated as argon in every set
            ExvInfo{.type=ff_t::OTHER, .Traube=volume(vdw::Ar), .Voronoi_implicit_H=volume(vdw::Ar), .MinimumFluctuation_implicit_H=volume(vdw::Ar), .Voronoi_explicit_H=volume(vdw::Ar), .MinimumFluctuation_explicit_H=volume(vdw::Ar), .vdw=volume(vdw::Ar)},
        };

        /**
         * @brief Extract a single volume set (column) from the descriptor table.
         */
        constexpr ExvSet make_set(std::optional<double> ExvInfo::* column) {
            ExvSet set{};
            for (const auto& row : exv_info_table) {
                if (set.volumes[static_cast<int>(row.type)].has_value()) {
                    throw ausaxs::except::runtime_error("constants::exv::detail::make_set: Duplicate row in exv_info_table.");
                }
                set.volumes[static_cast<int>(row.type)] = row.*column;
            }
            return set;
        }
    }

    constexpr detail::ExvSet Traube                        = detail::make_set(&detail::ExvInfo::Traube);
    constexpr detail::ExvSet Voronoi_implicit_H            = detail::make_set(&detail::ExvInfo::Voronoi_implicit_H);
    constexpr detail::ExvSet MinimumFluctuation_implicit_H = detail::make_set(&detail::ExvInfo::MinimumFluctuation_implicit_H);
    constexpr detail::ExvSet Voronoi_explicit_H            = detail::make_set(&detail::ExvInfo::Voronoi_explicit_H);
    constexpr detail::ExvSet MinimumFluctuation_explicit_H = detail::make_set(&detail::ExvInfo::MinimumFluctuation_explicit_H);
    constexpr detail::ExvSet vdw                           = detail::make_set(&detail::ExvInfo::vdw);

    // the fallback type and water must be present in every set
    static_assert(Traube.contains(ausaxs::form_factor::form_factor_t::OTHER) && Traube.contains(ausaxs::form_factor::form_factor_t::WATER));
    static_assert(Voronoi_implicit_H.contains(ausaxs::form_factor::form_factor_t::OTHER) && Voronoi_implicit_H.contains(ausaxs::form_factor::form_factor_t::WATER));
    static_assert(MinimumFluctuation_implicit_H.contains(ausaxs::form_factor::form_factor_t::OTHER) && MinimumFluctuation_implicit_H.contains(ausaxs::form_factor::form_factor_t::WATER));
    static_assert(Voronoi_explicit_H.contains(ausaxs::form_factor::form_factor_t::OTHER) && Voronoi_explicit_H.contains(ausaxs::form_factor::form_factor_t::WATER));
    static_assert(MinimumFluctuation_explicit_H.contains(ausaxs::form_factor::form_factor_t::OTHER) && MinimumFluctuation_explicit_H.contains(ausaxs::form_factor::form_factor_t::WATER));
    static_assert(vdw.contains(ausaxs::form_factor::form_factor_t::OTHER) && vdw.contains(ausaxs::form_factor::form_factor_t::WATER));

    constexpr double OH2 = 2.98*constexpr_math::pow(10., -23)*constexpr_math::pow(constants::SI::length::cm/constants::SI::length::A, 3);
    constexpr double Ar = detail::volume(constants::radius::vdw::Ar);
}
