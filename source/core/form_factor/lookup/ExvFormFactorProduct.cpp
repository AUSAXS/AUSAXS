// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <form_factor/lookup/ExvFormFactorProduct.h>
#include <form_factor/NormalizedFormFactor.h>
#include <settings/MoleculeSettings.h>
#include <constants/Constants.h>

#if CONSTEXPR_TABLES
    #define CONST constexpr
#else
    #define CONST const
#endif

using namespace ausaxs;
using namespace ausaxs::form_factor;

const FormFactorProduct& form_factor::lookup::exv::raw::get_product(unsigned int i, unsigned int j) {
    return get_table().index(i, j);
}

namespace {
    CONST form_factor::lookup::exv::table_t generate_exv_exv_table(const detail::ExvFormFactorSet& set) {
        form_factor::lookup::exv::table_t table;
        for (unsigned int i = 0; i < form_factor::get_count_without_excluded_volume(); ++i) {
            for (unsigned int j = 0; j < i; ++j) {
                table.index(i, j) = FormFactorProduct(
                    set.get(static_cast<form_factor_t>(i)), 
                    set.get(static_cast<form_factor_t>(j))
                );
                table.index(j, i) = table.index(i, j);
            }
            table.index(i, i) = FormFactorProduct(
                set.get(static_cast<form_factor_t>(i)), 
                set.get(static_cast<form_factor_t>(i))
            );
        }
        return table;
    }

    form_factor::lookup::exv::table_t custom_table_xx;
    auto ff_xx_default_table = generate_exv_exv_table(lookup::exv::standard);
}

const form_factor::lookup::exv::table_t& form_factor::lookup::exv::raw::get_table() {
    if (settings::molecule::exv_set == settings::molecule::ExvSet::Default) {
        return ff_xx_default_table;
    }
    switch (settings::molecule::exv_set) {
        case settings::molecule::ExvSet::Traube: {
            static auto table = generate_exv_exv_table(::detail::ExvFormFactorSet(constants::exv::Traube));
            return table;
        }
        case settings::molecule::ExvSet::Voronoi_explicit_H: {
            static auto table = generate_exv_exv_table(::detail::ExvFormFactorSet(constants::exv::Voronoi_explicit_H));
            return table;
        }
        case settings::molecule::ExvSet::Voronoi_implicit_H: {
            static auto table = generate_exv_exv_table(::detail::ExvFormFactorSet(constants::exv::Voronoi_implicit_H));
            return table;
        }
        case settings::molecule::ExvSet::MinimumFluctutation_explicit_H: {
            static auto table = generate_exv_exv_table(::detail::ExvFormFactorSet(constants::exv::MinimumFluctuation_explicit_H));
            return table;
        }
        case settings::molecule::ExvSet::MinimumFluctutation_implicit_H: {
            static auto table = generate_exv_exv_table(::detail::ExvFormFactorSet(constants::exv::MinimumFluctuation_implicit_H));
            return table;
        }
        case settings::molecule::ExvSet::vdw: {
            static auto table = generate_exv_exv_table(::detail::ExvFormFactorSet(constants::exv::vdw));
            return table;
        }
        case settings::molecule::ExvSet::Custom: {
            return custom_table_xx;
        }
    }
    return ff_xx_default_table;
}

const FormFactorProduct& form_factor::lookup::cross::raw::get_product(unsigned int i, unsigned int j) {
    return get_table().index(i, j);
}

namespace {
    CONST form_factor::lookup::cross::table_t generate_atom_exv_table(const detail::ExvFormFactorSet& set) {
        form_factor::lookup::cross::table_t table;
        for (unsigned int i = 0; i < form_factor::get_count_without_excluded_volume(); ++i) {
            for (unsigned int j = 0; j < form_factor::get_count_without_excluded_volume(); ++j) {
                table.index(i, j) = FormFactorProduct(
                    lookup::atomic::raw::get(static_cast<form_factor_t>(i)), 
                    set.get(static_cast<form_factor_t>(j))
                );
            }
        }
        return table;
    }

    form_factor::lookup::cross::table_t custom_table_ax;
    auto ff_ax_default_table = generate_atom_exv_table(lookup::exv::standard);
}

const form_factor::lookup::cross::table_t& form_factor::lookup::cross::raw::get_table() {
    if (settings::molecule::exv_set == settings::molecule::ExvSet::Default) {
        return ff_ax_default_table;
    }
    switch (settings::molecule::exv_set) {
        case settings::molecule::ExvSet::Traube: {
            static auto table = generate_atom_exv_table(::detail::ExvFormFactorSet(constants::exv::Traube));
            return table;
        }
        case settings::molecule::ExvSet::Voronoi_explicit_H: {
            static auto table = generate_atom_exv_table(::detail::ExvFormFactorSet(constants::exv::Voronoi_explicit_H));
            return table;
        }
        case settings::molecule::ExvSet::Voronoi_implicit_H: {
            static auto table = generate_atom_exv_table(::detail::ExvFormFactorSet(constants::exv::Voronoi_implicit_H));
            return table;
        }
        case settings::molecule::ExvSet::MinimumFluctutation_explicit_H: {
            static auto table = generate_atom_exv_table(::detail::ExvFormFactorSet(constants::exv::MinimumFluctuation_explicit_H));
            return table;
        }
        case settings::molecule::ExvSet::MinimumFluctutation_implicit_H: {
            static auto table = generate_atom_exv_table(::detail::ExvFormFactorSet(constants::exv::MinimumFluctuation_implicit_H));
            return table;
        }
        case settings::molecule::ExvSet::vdw: {
            static auto table = generate_atom_exv_table(::detail::ExvFormFactorSet(constants::exv::vdw));
            return table;
        }
        case settings::molecule::ExvSet::Custom: {
            return custom_table_ax;
        }
    }
    return ff_ax_default_table;
}

void form_factor::lookup::exv::raw::detail::set_custom_table(const constants::exv::detail::ExvSet& set) {
    ::detail::ExvFormFactorSet ff(set);
    custom_table_xx = generate_exv_exv_table(ff);
    custom_table_ax = generate_atom_exv_table(ff);
}