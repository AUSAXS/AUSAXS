// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <constants/Constants.h>
#include <form_factor/FormFactorTable.h>
#include <utility/Exceptions.h>

#include <array>
#include <string>
#include <string_view>

namespace ausaxs::form_factor {
    // The form factor type of an atom. This is intended to be used as an index for best performance.
    // NOLINTNEXTLINE(readability-enum-initial-value)
    enum class form_factor_t : int {
        EXCLUDED_VOLUME,    // excluded volume
        WATER,              // water
        OH = WATER,         // neutral oxygen with hydrogen
        H,                  // neutral hydrogen
        C,                  // neutral carbon
        CH,                 // neutral carbon with hydrogen
        CH2,                // neutral carbon with two hydrogens
        CH3,                // neutral carbon with three hydrogens
        N,                  // neutral nitrogen
        NH,                 // neutral nitrogen with hydrogen
        NH2,                // neutral nitrogen with two hydrogens
        NH3,                // neutral nitrogen with three hydrogens
        O,                  // neutral oxygen
        S,                  // neutral sulfur
        SH,                 // neutral sulfur with hydrogen
        OTHER,              // all other atoms
        COUNT,              // this will have the numerical value of the number of form factor types, and can thus be used to allocate arrays
        UNKNOWN,            // this is used to indicate that the form factor is unknown
    };

    /**
     * @brief The number of defined form factor types (including excluded volume).
     */
    constexpr int total_ff_count = static_cast<int>(form_factor::form_factor_t::COUNT);

    constexpr int exv_bin   = static_cast<int>(form_factor::form_factor_t::EXCLUDED_VOLUME);
    constexpr int water_bin = static_cast<int>(form_factor::form_factor_t::WATER);
    static_assert(exv_bin == 0, "form_factor::form_factor_t::EXCLUDED_VOLUME must be at index 0");
    static_assert(water_bin == 1, "form_factor::form_factor_t::WATER must be at index 1");

    constexpr int start_index_for_explicit_exv() {return static_cast<int>(form_factor::form_factor_t::EXCLUDED_VOLUME)+1;}

    namespace detail {
        /**
         * @brief The number of form factor slots currently in use.
         */
        inline int active_ff_count = total_ff_count;
    }

    /**
     * @brief Get the number of active form factors.
     */
    inline int get_active_count() noexcept {return detail::active_ff_count;}

    /**
     * @brief Descriptor of a single form factor type. 
     */
    struct FormFactorInfo {
        form_factor_t type;                                 // The type described. Must match the row index in the table below.
        std::string_view name;                              // Short name, used for I/O.
        constants::atom_t element;                          // The (heavy) atom of the type.
        int hydrogens;                                      // The number of implicit hydrogens bound to the element.
        int electrons;                                      // The number of electrons, including those of the bound hydrogens.
        double mass;                                        // The mass in amu, including the bound hydrogens.
        constants::form_factor::FiveGaussian coefficients;  // The five-Gaussian approximation of the vacuum form factor.
    };

    namespace detail {
        namespace ff = constants::form_factor;
        using constants::atom_t;
        using constants::mass::get_mass;

        /**
         * @brief The descriptor table of all form factor types, indexed by form_factor_t.
         *        To add a new type, add it to the enum and give it a row here. 
         *        Its excluded volume is described separately in form_factor/ExvTable.h, and is optional. 
         */
        constexpr std::array<FormFactorInfo, total_ff_count> ff_info_table = {{
            {.type=form_factor_t::EXCLUDED_VOLUME, .name="EXV", .element=atom_t::unknown, .hydrogens=0, .electrons=0,  .mass=0,                   .coefficients=ff::excluded_volume},
            {.type=form_factor_t::OH,              .name="OH",  .element=atom_t::O,       .hydrogens=1, .electrons=9,  .mass=16.999,              .coefficients=ff::OH_alc         },
            {.type=form_factor_t::H,               .name="H",   .element=atom_t::H,       .hydrogens=0, .electrons=1,  .mass=get_mass(atom_t::H), .coefficients=ff::H              },
            {.type=form_factor_t::C,               .name="C",   .element=atom_t::C,       .hydrogens=0, .electrons=6,  .mass=get_mass(atom_t::C), .coefficients=ff::C              },
            {.type=form_factor_t::CH,              .name="CH",  .element=atom_t::C,       .hydrogens=1, .electrons=7,  .mass=13.019,              .coefficients=ff::CH_sp3         },
            {.type=form_factor_t::CH2,             .name="CH2", .element=atom_t::C,       .hydrogens=2, .electrons=8,  .mass=14.027,              .coefficients=ff::CH2_sp3        },
            {.type=form_factor_t::CH3,             .name="CH3", .element=atom_t::C,       .hydrogens=3, .electrons=9,  .mass=15.035,              .coefficients=ff::CH3_sp3        },
            {.type=form_factor_t::N,               .name="N",   .element=atom_t::N,       .hydrogens=0, .electrons=7,  .mass=14.00674,            .coefficients=ff::N              },
            {.type=form_factor_t::NH,              .name="NH",  .element=atom_t::N,       .hydrogens=1, .electrons=8,  .mass=15.01474,            .coefficients=ff::NH             },
            {.type=form_factor_t::NH2,             .name="NH2", .element=atom_t::N,       .hydrogens=2, .electrons=9,  .mass=16.02274,            .coefficients=ff::NH2            },
            {.type=form_factor_t::NH3,             .name="NH3", .element=atom_t::N,       .hydrogens=3, .electrons=10, .mass=17.03074,            .coefficients=ff::NH3_plus       },
            {.type=form_factor_t::O,               .name="O",   .element=atom_t::O,       .hydrogens=0, .electrons=8,  .mass=15.999,              .coefficients=ff::O              },
            {.type=form_factor_t::S,               .name="S",   .element=atom_t::S,       .hydrogens=0, .electrons=16, .mass=32.06,               .coefficients=ff::S              },
            {.type=form_factor_t::SH,              .name="SH",  .element=atom_t::S,       .hydrogens=1, .electrons=17, .mass=33.06,               .coefficients=ff::SH             },
            {.type=form_factor_t::OTHER,           .name="OTH", .element=atom_t::Ar,      .hydrogens=0, .electrons=18, .mass=39.948,              .coefficients=ff::other          },
        }};

        constexpr bool ff_info_table_is_ordered() {
            for (int i = 0; i < total_ff_count; ++i) {
                if (static_cast<int>(ff_info_table[i].type) != i) {return false;}
            }
            return true;
        }
        static_assert(ff_info_table_is_ordered(), "form_factor::detail::ff_info_table must be ordered by form_factor_t.");

        constexpr bool is_tabulated(form_factor_t type) {
            return 0 <= static_cast<int>(type) && static_cast<int>(type) < total_ff_count;
        }

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
        if (type == ausaxs::form_factor::form_factor_t::COUNT) {return 0;}
        return ausaxs::form_factor::get_info(type).mass;
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
        return ausaxs::form_factor::get_info(type).electrons;
    }
}

namespace ausaxs::constants::charge {
    /**
     * @brief Get the effective charge based on the form factor evaluated at q=0.
     *        This represents the scattering power of the atom/group.
     */
    double get_ff_charge(ausaxs::form_factor::form_factor_t type);

    /**
     * @brief Get the effective charge of an atom, with its element as a fallback.
     *        When the form factor is unknown, the atomic charge is returned instead.
     */
    double get_ff_charge(ausaxs::form_factor::form_factor_t type, ausaxs::constants::atom_t fallback_element);
}