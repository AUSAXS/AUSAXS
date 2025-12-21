// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <form_factor/lookup/ExvFormFactorProduct.h>
#include <form_factor/lookup/FormFactorProductBase.h>
#include <form_factor/FormFactor.h>
#include <settings/MoleculeSettings.h>
#include <constants/Constants.h>

using namespace ausaxs;
using namespace ausaxs::form_factor;

namespace {
    struct RawFormFactorLookup {
        static constexpr const FormFactor& get(form_factor_t type) {
            return lookup::atomic::raw::get(type);
        }
    };

    form_factor::lookup::exv::table_t custom_table_xx;
    form_factor::lookup::cross::table_t custom_table_ax;
    auto ff_xx_default_table = lookup::detail::generate_exv_table<lookup::exv::table_t>(lookup::exv::standard);
    auto ff_ax_default_table = lookup::detail::generate_cross_table<lookup::cross::table_t, RawFormFactorLookup>(lookup::exv::standard);
}

const FormFactorProduct& form_factor::lookup::exv::raw::get_product(unsigned int i, unsigned int j) {
    return get_table().index(i, j);
}

const form_factor::lookup::exv::table_t& form_factor::lookup::exv::raw::get_table() {
    return lookup::detail::get_exv_table_for_settings(
        ff_xx_default_table, custom_table_xx,
        [](const form_factor::detail::ExvFormFactorSet& set) { return lookup::detail::generate_exv_table<lookup::exv::table_t>(set); }
    );
}

const FormFactorProduct& form_factor::lookup::cross::raw::get_product(unsigned int i, unsigned int j) {
    return get_table().index(i, j);
}

const form_factor::lookup::cross::table_t& form_factor::lookup::cross::raw::get_table() {
    return lookup::detail::get_exv_table_for_settings(
        ff_ax_default_table, custom_table_ax,
        [](const form_factor::detail::ExvFormFactorSet& set) { return lookup::detail::generate_cross_table<lookup::cross::table_t, RawFormFactorLookup>(set); }
    );
}

void form_factor::lookup::exv::raw::detail::set_custom_table(const constants::exv::detail::ExvSet& set) {
    form_factor::detail::ExvFormFactorSet ff(set);
    custom_table_xx = lookup::detail::generate_exv_table<lookup::exv::table_t>(ff);
    custom_table_ax = lookup::detail::generate_cross_table<lookup::cross::table_t, RawFormFactorLookup>(ff);
}