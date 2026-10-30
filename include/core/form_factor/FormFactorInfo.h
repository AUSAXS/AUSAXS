// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <constants/Constants.h>
#include <form_factor/FormFactorType.h>
#include <form_factor/GroupGeometry.h>
#include <form_factor/neutron/NeutronFormFactorTable.h>
#include <form_factor/xray/XrayFormFactorTable.h>
#include <settings/ScatteringSettings.h>
#include <utility/Exceptions.h>

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
        group::Geometry geometry;                           // The implicit hydrogens bound to the element.
        xray::coefficients::FiveGaussian xray_coefficients; // The five-Gaussian approximation of the vacuum X-ray form factor of the element.
        double neutron_scattering_length;                   // The coherent neutron scattering length of the element in fm.
    };

    namespace detail {
        namespace ff = xray::coefficients;
        namespace b = neutron::scattering_length;
        using constants::atom_t;

        /**
         * @brief The descriptor table of all form factor types, indexed by form_factor_t.
         *        To add a new type, add it to the enum and give it a row here. 
         *        Its excluded volume is described separately in form_factor/ExvTable.h, and is optional. 
         */
        constexpr std::array<FormFactorInfo, total_ff_count> ff_info_table = {{
            {.type=form_factor_t::EXCLUDED_VOLUME, .name="EXV", .element=atom_t::unknown, .geometry=group::atom, .xray_coefficients=ff::excluded_volume, .neutron_scattering_length=b::excluded_volume},
            {.type=form_factor_t::OH,              .name="OH",  .element=atom_t::O,       .geometry=group::OH,   .xray_coefficients=ff::O,               .neutron_scattering_length=b::O              },
            {.type=form_factor_t::C,               .name="C",   .element=atom_t::C,       .geometry=group::atom, .xray_coefficients=ff::C,               .neutron_scattering_length=b::C              },
            {.type=form_factor_t::CH,              .name="CH",  .element=atom_t::C,       .geometry=group::CH,   .xray_coefficients=ff::C,               .neutron_scattering_length=b::C              },
            {.type=form_factor_t::CH2,             .name="CH2", .element=atom_t::C,       .geometry=group::CH2,  .xray_coefficients=ff::C,               .neutron_scattering_length=b::C              },
            {.type=form_factor_t::CH3,             .name="CH3", .element=atom_t::C,       .geometry=group::CH3,  .xray_coefficients=ff::C,               .neutron_scattering_length=b::C              },
            {.type=form_factor_t::N,               .name="N",   .element=atom_t::N,       .geometry=group::atom, .xray_coefficients=ff::N,               .neutron_scattering_length=b::N              },
            {.type=form_factor_t::NH,              .name="NH",  .element=atom_t::N,       .geometry=group::NH,   .xray_coefficients=ff::N,               .neutron_scattering_length=b::N              },
            {.type=form_factor_t::NH2,             .name="NH2", .element=atom_t::N,       .geometry=group::NH2,  .xray_coefficients=ff::N,               .neutron_scattering_length=b::N              },
            {.type=form_factor_t::NH3,             .name="NH3", .element=atom_t::N,       .geometry=group::NH3,  .xray_coefficients=ff::N,               .neutron_scattering_length=b::N              },
            {.type=form_factor_t::O,               .name="O",   .element=atom_t::O,       .geometry=group::atom, .xray_coefficients=ff::O,               .neutron_scattering_length=b::O              },
            {.type=form_factor_t::S,               .name="S",   .element=atom_t::S,       .geometry=group::atom, .xray_coefficients=ff::S,               .neutron_scattering_length=b::S              },
            {.type=form_factor_t::SH,              .name="SH",  .element=atom_t::S,       .geometry=group::SH,   .xray_coefficients=ff::S,               .neutron_scattering_length=b::S              },
            {.type=form_factor_t::H,               .name="H",   .element=atom_t::H,       .geometry=group::atom, .xray_coefficients=ff::H,               .neutron_scattering_length=b::H              },
            {.type=form_factor_t::D,               .name="D",   .element=atom_t::D,       .geometry=group::atom, .xray_coefficients=ff::H,               .neutron_scattering_length=b::D              },
            {.type=form_factor_t::He,              .name="He",  .element=atom_t::He,      .geometry=group::atom, .xray_coefficients=ff::He,              .neutron_scattering_length=b::He             },
            {.type=form_factor_t::Li,              .name="Li",  .element=atom_t::Li,      .geometry=group::atom, .xray_coefficients=ff::Li,              .neutron_scattering_length=b::Li             },
            {.type=form_factor_t::Be,              .name="Be",  .element=atom_t::Be,      .geometry=group::atom, .xray_coefficients=ff::Be,              .neutron_scattering_length=b::Be             },
            {.type=form_factor_t::B,               .name="B",   .element=atom_t::B,       .geometry=group::atom, .xray_coefficients=ff::B,               .neutron_scattering_length=b::B              },
            {.type=form_factor_t::F,               .name="F",   .element=atom_t::F,       .geometry=group::atom, .xray_coefficients=ff::F,               .neutron_scattering_length=b::F              },
            {.type=form_factor_t::Ne,              .name="Ne",  .element=atom_t::Ne,      .geometry=group::atom, .xray_coefficients=ff::Ne,              .neutron_scattering_length=b::Ne             },
            {.type=form_factor_t::Na,              .name="Na",  .element=atom_t::Na,      .geometry=group::atom, .xray_coefficients=ff::Na,              .neutron_scattering_length=b::Na             },
            {.type=form_factor_t::Mg,              .name="Mg",  .element=atom_t::Mg,      .geometry=group::atom, .xray_coefficients=ff::Mg,              .neutron_scattering_length=b::Mg             },
            {.type=form_factor_t::Al,              .name="Al",  .element=atom_t::Al,      .geometry=group::atom, .xray_coefficients=ff::Al,              .neutron_scattering_length=b::Al             },
            {.type=form_factor_t::Si,              .name="Si",  .element=atom_t::Si,      .geometry=group::atom, .xray_coefficients=ff::Si,              .neutron_scattering_length=b::Si             },
            {.type=form_factor_t::P,               .name="P",   .element=atom_t::P,       .geometry=group::atom, .xray_coefficients=ff::P,               .neutron_scattering_length=b::P              },
            {.type=form_factor_t::Cl,              .name="Cl",  .element=atom_t::Cl,      .geometry=group::atom, .xray_coefficients=ff::Cl,              .neutron_scattering_length=b::Cl             },
            {.type=form_factor_t::Ar,              .name="Ar",  .element=atom_t::Ar,      .geometry=group::atom, .xray_coefficients=ff::Ar,              .neutron_scattering_length=b::Ar             },
            {.type=form_factor_t::K,               .name="K",   .element=atom_t::K,       .geometry=group::atom, .xray_coefficients=ff::K,               .neutron_scattering_length=b::K              },
            {.type=form_factor_t::Ca,              .name="Ca",  .element=atom_t::Ca,      .geometry=group::atom, .xray_coefficients=ff::Ca,              .neutron_scattering_length=b::Ca             },
            {.type=form_factor_t::Sc,              .name="Sc",  .element=atom_t::Sc,      .geometry=group::atom, .xray_coefficients=ff::Sc,              .neutron_scattering_length=b::Sc             },
            {.type=form_factor_t::Ti,              .name="Ti",  .element=atom_t::Ti,      .geometry=group::atom, .xray_coefficients=ff::Ti,              .neutron_scattering_length=b::Ti             },
            {.type=form_factor_t::V,               .name="V",   .element=atom_t::V,       .geometry=group::atom, .xray_coefficients=ff::V,               .neutron_scattering_length=b::V              },
            {.type=form_factor_t::Cr,              .name="Cr",  .element=atom_t::Cr,      .geometry=group::atom, .xray_coefficients=ff::Cr,              .neutron_scattering_length=b::Cr             },
            {.type=form_factor_t::Mn,              .name="Mn",  .element=atom_t::Mn,      .geometry=group::atom, .xray_coefficients=ff::Mn,              .neutron_scattering_length=b::Mn             },
            {.type=form_factor_t::Fe,              .name="Fe",  .element=atom_t::Fe,      .geometry=group::atom, .xray_coefficients=ff::Fe,              .neutron_scattering_length=b::Fe             },
            {.type=form_factor_t::Co,              .name="Co",  .element=atom_t::Co,      .geometry=group::atom, .xray_coefficients=ff::Co,              .neutron_scattering_length=b::Co             },
            {.type=form_factor_t::Ni,              .name="Ni",  .element=atom_t::Ni,      .geometry=group::atom, .xray_coefficients=ff::Ni,              .neutron_scattering_length=b::Ni             },
            {.type=form_factor_t::Cu,              .name="Cu",  .element=atom_t::Cu,      .geometry=group::atom, .xray_coefficients=ff::Cu,              .neutron_scattering_length=b::Cu             },
            {.type=form_factor_t::Zn,              .name="Zn",  .element=atom_t::Zn,      .geometry=group::atom, .xray_coefficients=ff::Zn,              .neutron_scattering_length=b::Zn             },
            {.type=form_factor_t::Ga,              .name="Ga",  .element=atom_t::Ga,      .geometry=group::atom, .xray_coefficients=ff::Ga,              .neutron_scattering_length=b::Ga             },
            {.type=form_factor_t::Ge,              .name="Ge",  .element=atom_t::Ge,      .geometry=group::atom, .xray_coefficients=ff::Ge,              .neutron_scattering_length=b::Ge             },
            {.type=form_factor_t::As,              .name="As",  .element=atom_t::As,      .geometry=group::atom, .xray_coefficients=ff::As,              .neutron_scattering_length=b::As             },
            {.type=form_factor_t::Se,              .name="Se",  .element=atom_t::Se,      .geometry=group::atom, .xray_coefficients=ff::Se,              .neutron_scattering_length=b::Se             },
            {.type=form_factor_t::Br,              .name="Br",  .element=atom_t::Br,      .geometry=group::atom, .xray_coefficients=ff::Br,              .neutron_scattering_length=b::Br             },
            {.type=form_factor_t::Kr,              .name="Kr",  .element=atom_t::Kr,      .geometry=group::atom, .xray_coefficients=ff::Kr,              .neutron_scattering_length=b::Kr             },
            {.type=form_factor_t::Rb,              .name="Rb",  .element=atom_t::Rb,      .geometry=group::atom, .xray_coefficients=ff::Rb,              .neutron_scattering_length=b::Rb             },
            {.type=form_factor_t::Sr,              .name="Sr",  .element=atom_t::Sr,      .geometry=group::atom, .xray_coefficients=ff::Sr,              .neutron_scattering_length=b::Sr             },
            {.type=form_factor_t::Y,               .name="Y",   .element=atom_t::Y,       .geometry=group::atom, .xray_coefficients=ff::Y,               .neutron_scattering_length=b::Y              },
            {.type=form_factor_t::Zr,              .name="Zr",  .element=atom_t::Zr,      .geometry=group::atom, .xray_coefficients=ff::Zr,              .neutron_scattering_length=b::Zr             },
            {.type=form_factor_t::Nb,              .name="Nb",  .element=atom_t::Nb,      .geometry=group::atom, .xray_coefficients=ff::Nb,              .neutron_scattering_length=b::Nb             },
            {.type=form_factor_t::Mo,              .name="Mo",  .element=atom_t::Mo,      .geometry=group::atom, .xray_coefficients=ff::Mo,              .neutron_scattering_length=b::Mo             },
            {.type=form_factor_t::Tc,              .name="Tc",  .element=atom_t::Tc,      .geometry=group::atom, .xray_coefficients=ff::Tc,              .neutron_scattering_length=b::Tc             },
            {.type=form_factor_t::Ru,              .name="Ru",  .element=atom_t::Ru,      .geometry=group::atom, .xray_coefficients=ff::Ru,              .neutron_scattering_length=b::Ru             },
            {.type=form_factor_t::Rh,              .name="Rh",  .element=atom_t::Rh,      .geometry=group::atom, .xray_coefficients=ff::Rh,              .neutron_scattering_length=b::Rh             },
            {.type=form_factor_t::Pd,              .name="Pd",  .element=atom_t::Pd,      .geometry=group::atom, .xray_coefficients=ff::Pd,              .neutron_scattering_length=b::Pd             },
            {.type=form_factor_t::Ag,              .name="Ag",  .element=atom_t::Ag,      .geometry=group::atom, .xray_coefficients=ff::Ag,              .neutron_scattering_length=b::Ag             },
            {.type=form_factor_t::Cd,              .name="Cd",  .element=atom_t::Cd,      .geometry=group::atom, .xray_coefficients=ff::Cd,              .neutron_scattering_length=b::Cd             },
            {.type=form_factor_t::In,              .name="In",  .element=atom_t::In,      .geometry=group::atom, .xray_coefficients=ff::In,              .neutron_scattering_length=b::In             },
            {.type=form_factor_t::Sn,              .name="Sn",  .element=atom_t::Sn,      .geometry=group::atom, .xray_coefficients=ff::Sn,              .neutron_scattering_length=b::Sn             },
            {.type=form_factor_t::Sb,              .name="Sb",  .element=atom_t::Sb,      .geometry=group::atom, .xray_coefficients=ff::Sb,              .neutron_scattering_length=b::Sb             },
            {.type=form_factor_t::Te,              .name="Te",  .element=atom_t::Te,      .geometry=group::atom, .xray_coefficients=ff::Te,              .neutron_scattering_length=b::Te             },
            {.type=form_factor_t::I,               .name="I",   .element=atom_t::I,       .geometry=group::atom, .xray_coefficients=ff::I,               .neutron_scattering_length=b::I              },
            {.type=form_factor_t::Xe,              .name="Xe",  .element=atom_t::Xe,      .geometry=group::atom, .xray_coefficients=ff::Xe,              .neutron_scattering_length=b::Xe             },
            {.type=form_factor_t::Cs,              .name="Cs",  .element=atom_t::Cs,      .geometry=group::atom, .xray_coefficients=ff::Cs,              .neutron_scattering_length=b::Cs             },
            {.type=form_factor_t::Ba,              .name="Ba",  .element=atom_t::Ba,      .geometry=group::atom, .xray_coefficients=ff::Ba,              .neutron_scattering_length=b::Ba             },
            {.type=form_factor_t::Hf,              .name="Hf",  .element=atom_t::Hf,      .geometry=group::atom, .xray_coefficients=ff::Hf,              .neutron_scattering_length=b::Hf             },
            {.type=form_factor_t::Ta,              .name="Ta",  .element=atom_t::Ta,      .geometry=group::atom, .xray_coefficients=ff::Ta,              .neutron_scattering_length=b::Ta             },
            {.type=form_factor_t::W,               .name="W",   .element=atom_t::W,       .geometry=group::atom, .xray_coefficients=ff::W,               .neutron_scattering_length=b::W              },
            {.type=form_factor_t::Re,              .name="Re",  .element=atom_t::Re,      .geometry=group::atom, .xray_coefficients=ff::Re,              .neutron_scattering_length=b::Re             },
            {.type=form_factor_t::Os,              .name="Os",  .element=atom_t::Os,      .geometry=group::atom, .xray_coefficients=ff::Os,              .neutron_scattering_length=b::Os             },
            {.type=form_factor_t::Ir,              .name="Ir",  .element=atom_t::Ir,      .geometry=group::atom, .xray_coefficients=ff::Ir,              .neutron_scattering_length=b::Ir             },
            {.type=form_factor_t::Pt,              .name="Pt",  .element=atom_t::Pt,      .geometry=group::atom, .xray_coefficients=ff::Pt,              .neutron_scattering_length=b::Pt             },
            {.type=form_factor_t::Au,              .name="Au",  .element=atom_t::Au,      .geometry=group::atom, .xray_coefficients=ff::Au,              .neutron_scattering_length=b::Au             },
            {.type=form_factor_t::Tl,              .name="Tl",  .element=atom_t::Tl,      .geometry=group::atom, .xray_coefficients=ff::Tl,              .neutron_scattering_length=b::Tl             },
            {.type=form_factor_t::Pb,              .name="Pb",  .element=atom_t::Pb,      .geometry=group::atom, .xray_coefficients=ff::Pb,              .neutron_scattering_length=b::Pb             },
            {.type=form_factor_t::Bi,              .name="Bi",  .element=atom_t::Bi,      .geometry=group::atom, .xray_coefficients=ff::Bi,              .neutron_scattering_length=b::Bi             },
            {.type=form_factor_t::OTHER,           .name="OTH", .element=atom_t::Ar,      .geometry=group::atom, .xray_coefficients=ff::other,           .neutron_scattering_length=b::other          },
        }};

        constexpr bool ff_info_table_is_ordered() {
            for (int i = 0; i < total_ff_count; ++i) {
                if (static_cast<int>(ff_info_table[i].type) != i) {return false;}
            }
            return true;
        }
        static_assert(ff_info_table_is_ordered(), "form_factor::detail::ff_info_table must be ordered by form_factor_t.");

        /**
         * @brief Map from an element to the form factor type of the bare element. Elements without their own type map to OTHER.
         */
        constexpr auto element_table = [] () {
            std::array<form_factor_t, static_cast<int>(atom_t::unknown)+1> table;
            table.fill(form_factor_t::OTHER);
            for (int i = start_index_for_explicit_exv(); i < total_ff_count; ++i) {
                const auto& info = ff_info_table[i];
                if (info.geometry.hydrogens == 0 && info.element != atom_t::unknown && info.type != form_factor_t::OTHER) {table[static_cast<int>(info.element)] = info.type;}
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
        return get_info(type).geometry.hydrogens != 0;
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
        return get_mass(info.element) + info.geometry.hydrogens*get_mass(atom_t::H);
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
        return get_charge(info.element) + info.geometry.hydrogens;
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
