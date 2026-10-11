// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <constants/Constants.h>
#include <form_factor/FormFactorType.h>
#include <form_factor/neutron/NeutronFormFactorTable.h>
#include <form_factor/xray/XrayFormFactorTable.h>
#include <settings/ScatteringSettings.h>
#include <utility/Exceptions.h>

#include <algorithm>
#include <array>
#include <string>
#include <string_view>

namespace ausaxs::form_factor {
    /**
     * @brief Descriptor of a single form factor type. 
     */
    struct FormFactorInfo {
        form_factor_t type;                                 // The type described. Must match the row index in the table below.
        std::string_view name;                              // Short name, used for I/O.
        constants::atom_t element;                          // The (heavy) atom of the type.
        int hydrogens;                                      // The number of implicit hydrogens bound to the element.
        xray::coefficients::FiveGaussian xray_coefficients; // The five-Gaussian approximation of the vacuum X-ray form factor.
        neutron::coefficients::Group neutron_coefficients;  // The scattering length and hydrogen geometry of the neutron form factor.
    };

    namespace detail {
        namespace ff = xray::coefficients;
        namespace nff = neutron::coefficients;
        using constants::atom_t;

        /**
         * @brief The descriptor table of all form factor types, indexed by form_factor_t.
         *        To add a new type, add it to the enum and give it a row here. 
         *        Its excluded volume is described separately in form_factor/ExvTable.h, and is optional. 
         */
        constexpr std::array<FormFactorInfo, total_ff_count> ff_info_table = {{
            {.type=form_factor_t::EXCLUDED_VOLUME, .name="EXV", .element=atom_t::unknown, .hydrogens=0, .xray_coefficients=ff::excluded_volume, .neutron_coefficients=nff::excluded_volume},
            {.type=form_factor_t::OH,              .name="OH",  .element=atom_t::O,       .hydrogens=1, .xray_coefficients=ff::OH_alc,          .neutron_coefficients=nff::OH             },
            {.type=form_factor_t::H,               .name="H",   .element=atom_t::H,       .hydrogens=0, .xray_coefficients=ff::H,               .neutron_coefficients=nff::H              },
            {.type=form_factor_t::C,               .name="C",   .element=atom_t::C,       .hydrogens=0, .xray_coefficients=ff::C,               .neutron_coefficients=nff::C              },
            {.type=form_factor_t::CH,              .name="CH",  .element=atom_t::C,       .hydrogens=1, .xray_coefficients=ff::CH_sp3,          .neutron_coefficients=nff::CH             },
            {.type=form_factor_t::CH2,             .name="CH2", .element=atom_t::C,       .hydrogens=2, .xray_coefficients=ff::CH2_sp3,         .neutron_coefficients=nff::CH2            },
            {.type=form_factor_t::CH3,             .name="CH3", .element=atom_t::C,       .hydrogens=3, .xray_coefficients=ff::CH3_sp3,         .neutron_coefficients=nff::CH3            },
            {.type=form_factor_t::N,               .name="N",   .element=atom_t::N,       .hydrogens=0, .xray_coefficients=ff::N,               .neutron_coefficients=nff::N              },
            {.type=form_factor_t::NH,              .name="NH",  .element=atom_t::N,       .hydrogens=1, .xray_coefficients=ff::NH,              .neutron_coefficients=nff::NH             },
            {.type=form_factor_t::NH2,             .name="NH2", .element=atom_t::N,       .hydrogens=2, .xray_coefficients=ff::NH2,             .neutron_coefficients=nff::NH2            },
            {.type=form_factor_t::NH3,             .name="NH3", .element=atom_t::N,       .hydrogens=3, .xray_coefficients=ff::NH3_plus,        .neutron_coefficients=nff::NH3            },
            {.type=form_factor_t::O,               .name="O",   .element=atom_t::O,       .hydrogens=0, .xray_coefficients=ff::O,               .neutron_coefficients=nff::O              },
            {.type=form_factor_t::S,               .name="S",   .element=atom_t::S,       .hydrogens=0, .xray_coefficients=ff::S,               .neutron_coefficients=nff::S              },
            {.type=form_factor_t::SH,              .name="SH",  .element=atom_t::S,       .hydrogens=1, .xray_coefficients=ff::SH,              .neutron_coefficients=nff::SH             },
            {.type=form_factor_t::F,               .name="F",   .element=atom_t::F,       .hydrogens=0, .xray_coefficients=ff::F,               .neutron_coefficients=nff::F              },
            {.type=form_factor_t::Na,              .name="Na",  .element=atom_t::Na,      .hydrogens=0, .xray_coefficients=ff::Na,              .neutron_coefficients=nff::Na             },
            {.type=form_factor_t::Mg,              .name="Mg",  .element=atom_t::Mg,      .hydrogens=0, .xray_coefficients=ff::Mg,              .neutron_coefficients=nff::Mg             },
            {.type=form_factor_t::P,               .name="P",   .element=atom_t::P,       .hydrogens=0, .xray_coefficients=ff::P,               .neutron_coefficients=nff::P              },
            {.type=form_factor_t::Cl,              .name="Cl",  .element=atom_t::Cl,      .hydrogens=0, .xray_coefficients=ff::Cl,              .neutron_coefficients=nff::Cl             },
            {.type=form_factor_t::K,               .name="K",   .element=atom_t::K,       .hydrogens=0, .xray_coefficients=ff::K,               .neutron_coefficients=nff::K              },
            {.type=form_factor_t::Ca,              .name="Ca",  .element=atom_t::Ca,      .hydrogens=0, .xray_coefficients=ff::Ca,              .neutron_coefficients=nff::Ca             },
            {.type=form_factor_t::Mn,              .name="Mn",  .element=atom_t::Mn,      .hydrogens=0, .xray_coefficients=ff::Mn,              .neutron_coefficients=nff::Mn             },
            {.type=form_factor_t::Fe,              .name="Fe",  .element=atom_t::Fe,      .hydrogens=0, .xray_coefficients=ff::Fe,              .neutron_coefficients=nff::Fe             },
            {.type=form_factor_t::Co,              .name="Co",  .element=atom_t::Co,      .hydrogens=0, .xray_coefficients=ff::Co,              .neutron_coefficients=nff::Co             },
            {.type=form_factor_t::Ni,              .name="Ni",  .element=atom_t::Ni,      .hydrogens=0, .xray_coefficients=ff::Ni,              .neutron_coefficients=nff::Ni             },
            {.type=form_factor_t::Cu,              .name="Cu",  .element=atom_t::Cu,      .hydrogens=0, .xray_coefficients=ff::Cu,              .neutron_coefficients=nff::Cu             },
            {.type=form_factor_t::Zn,              .name="Zn",  .element=atom_t::Zn,      .hydrogens=0, .xray_coefficients=ff::Zn,              .neutron_coefficients=nff::Zn             },
            {.type=form_factor_t::Se,              .name="Se",  .element=atom_t::Se,      .hydrogens=0, .xray_coefficients=ff::Se,              .neutron_coefficients=nff::Se             },
            {.type=form_factor_t::Br,              .name="Br",  .element=atom_t::Br,      .hydrogens=0, .xray_coefficients=ff::Br,              .neutron_coefficients=nff::Br             },
            {.type=form_factor_t::I,               .name="I",   .element=atom_t::I,       .hydrogens=0, .xray_coefficients=ff::I,               .neutron_coefficients=nff::I              },
            {.type=form_factor_t::OTHER,           .name="OTH", .element=atom_t::Ar,      .hydrogens=0, .xray_coefficients=ff::other,           .neutron_coefficients=nff::other          },
        }};

        constexpr bool ff_info_table_is_ordered() {
            for (int i = 0; i < total_ff_count; ++i) {
                if (static_cast<int>(ff_info_table[i].type) != i) {return false;}
            }
            return true;
        }
        static_assert(ff_info_table_is_ordered(), "form_factor::detail::ff_info_table must be ordered by form_factor_t.");

        constexpr bool ff_info_table_hydrogens_agree() {
            return std::ranges::all_of(ff_info_table, [](const auto& info) {
                return info.neutron_coefficients.hydrogens == info.hydrogens;
            });
        }
        static_assert(ff_info_table_hydrogens_agree(), "form_factor::detail::ff_info_table: the neutron coefficients of each row must have the same number of hydrogens as the row.");

        /**
         * @brief Map from an element to the form factor type of the bare element. Elements without their own type map to OTHER.
         */
        constexpr auto element_table = [] () {
            std::array<form_factor_t, static_cast<int>(atom_t::unknown)+1> table;
            table.fill(form_factor_t::OTHER);
            for (int i = start_index_for_explicit_exv(); i < total_ff_count; ++i) {
                const auto& info = ff_info_table[i];
                if (info.hydrogens == 0 && info.element != atom_t::unknown) {table[static_cast<int>(info.element)] = info.type;}
            }
            return table;
        }();
    }

    /**
     * @brief Get the descriptor of a form factor type.
     */
    constexpr const FormFactorInfo& get_info(form_factor_t type) {
        if (!detail::is_tabulated(type)) {
            throw ausaxs::except::runtime_error("form_factor::get_info: Invalid form factor type (enum " + std::to_string(static_cast<int>(type)) + ")");
        }
        return detail::ff_info_table[static_cast<int>(type)];
    }

    /**
     * @brief Check if a form factor type represents an atomic group, i.e. already includes the hydrogens bound to the element.
     */
    constexpr bool has_implicit_hydrogens(form_factor_t type) {
        return get_info(type).hydrogens != 0;
    }

    [[maybe_unused]] static std::string to_string(form_factor_t type) {
        switch (type) {
            case form_factor_t::COUNT: return "CNT";
            case form_factor_t::UNKNOWN: return "UNK";
            default: return std::string(get_info(type).name);
        }
    }

    [[maybe_unused]] static form_factor_t from_string(const std::string& str) {
        for (const auto& info : detail::ff_info_table) {
            if (info.name == str) {return info.type;}
        }
        if (str == "CNT") return form_factor_t::COUNT;
        if (str == "UNK") return form_factor_t::UNKNOWN;
        throw ausaxs::except::runtime_error("form_factor::from_string: Unknown form factor string \"" + str + "\"");
    }

    /**
     * @brief Get the form factor type based on an atom type.
     *        In case the atom type is not recognized, the default form factor (argon) is returned.
     */
    constexpr form_factor_t get_type(constants::atom_t atom_type) {
        int i = static_cast<int>(atom_type);
        if (i < 0 || static_cast<int>(detail::element_table.size()) <= i) {return form_factor_t::OTHER;}
        return detail::element_table[i];
    }

    /**
     * @brief Get the form factor type based on an atom type and an atomic group.
     *        The atomic group takes priority. Only if the atomic group is not recognized, the atom type is used.
     *        In case either the atomic group or the atom type is not recognized, the default form factor (argon) is returned.
     */
    constexpr form_factor_t get_type(constants::atom_t atom_type, constants::atomic_group_t atomic_group) {
        switch(atomic_group) {
            case constants::atomic_group_t::CH: return form_factor_t::CH;
            case constants::atomic_group_t::CH2: return form_factor_t::CH2;
            case constants::atomic_group_t::CH3: return form_factor_t::CH3;
            case constants::atomic_group_t::NH: return form_factor_t::NH;
            case constants::atomic_group_t::NH2: return form_factor_t::NH2;
            case constants::atomic_group_t::NH3: return form_factor_t::NH3;
            case constants::atomic_group_t::OH: return form_factor_t::OH;
            case constants::atomic_group_t::SH: return form_factor_t::SH;
            default: return get_type(atom_type);
        }
    }

    constexpr constants::atom_t to_atom_type(form_factor_t ff_type) {
        if (!detail::is_tabulated(ff_type) || ff_type == form_factor_t::EXCLUDED_VOLUME) {return constants::atom_t::Ar;}
        return get_info(ff_type).element;
    }
}

namespace ausaxs::constants::mass {
    /**
    * @brief Get the mass of an atom in amu.
    */
    constexpr double get_mass(ausaxs::form_factor::form_factor_t type) {
        using ausaxs::form_factor::form_factor_t;
        if (type == form_factor_t::COUNT || type == form_factor_t::EXCLUDED_VOLUME) {return 0;}
        const auto& info = ausaxs::form_factor::get_info(type);
        return get_mass(info.element) + info.hydrogens*get_mass(atom_t::H);
    }
}

namespace ausaxs::constants::radius {
    inline double get_vdw_radius(ausaxs::form_factor::form_factor_t type) {
        using ausaxs::form_factor::form_factor_t;
        if (type == form_factor_t::UNKNOWN) {return get_vdw_radius(constants::atom_t::C);} // ensure unknown atoms still occupy space in the grid
        if (!ausaxs::form_factor::detail::is_tabulated(type) || type == form_factor_t::EXCLUDED_VOLUME) {
            throw ausaxs::except::runtime_error("constants::radius::get_vdw_radius: Unknown form factor type \"" + ausaxs::form_factor::to_string(type) + "\"");
        }
        return get_vdw_radius(ausaxs::form_factor::get_info(type).element);
    }
}

namespace ausaxs::constants::charge::nuclear {
    /**
     * @brief Get the charge of an atom in e.
     */
    constexpr int get_charge(ausaxs::form_factor::form_factor_t type) {
        if (type == ausaxs::form_factor::form_factor_t::EXCLUDED_VOLUME) {return 0;}
        const auto& info = ausaxs::form_factor::get_info(type);
        return get_charge(info.element) + info.hydrogens;
    }
}

namespace ausaxs::constants::charge {
    /**
     * @brief Get the effective charge based on the form factor evaluated at q=0.
     *        This represents the scattering power of the atom/group for the given radiation.
     */
    template<settings::scattering::Radiation radiation>
    double get_ff_charge(ausaxs::form_factor::form_factor_t type);

    /**
     * @brief Get the effective charge of an atom, with its element as a fallback.
     *        When the form factor is unknown, the atomic charge is returned instead.
     */
    template<settings::scattering::Radiation radiation>
    double get_ff_charge(ausaxs::form_factor::form_factor_t type, ausaxs::constants::atom_t fallback_element);

    /**
     * @brief Get the effective charge for the radiation of the current settings, with its element as a fallback.
     */
    double get_ff_charge(ausaxs::form_factor::form_factor_t type, ausaxs::constants::atom_t fallback_element);

    /**
     * @brief Get the effective charge for the radiation of the current settings.
     */
    double get_ff_charge(ausaxs::form_factor::form_factor_t type);
}
