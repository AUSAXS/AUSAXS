// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <form_factor/lookup/ExvFormFactorProduct.h>
#include <form_factor/lookup/CustomExvTable.h>
#include <form_factor/lookup/detail/FormFactorProductBase.h>
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

    auto ff_xx_default_table = lookup::detail::generate_exv_table<lookup::exv::table_t>(lookup::exv::standard);
    auto ff_ax_default_table = lookup::detail::generate_cross_table<lookup::cross::table_t, RawFormFactorLookup>(lookup::exv::standard);
}

const FormFactorProduct& form_factor::lookup::exv::raw::get_product(unsigned int i, unsigned int j) {
    return get_table().index(i, j);
}

const form_factor::lookup::exv::table_t& form_factor::lookup::exv::raw::get_table() {
    return lookup::detail::get_exv_table_for_settings(
        ff_xx_default_table, lookup::exv::detail::custom_table_xx_raw,
        [](const form_factor::detail::ExvFormFactorSet& set) { return lookup::detail::generate_exv_table<lookup::exv::table_t>(set); }
    );
}

const FormFactorProduct& form_factor::lookup::cross::raw::get_product(unsigned int i, unsigned int j) {
    return get_table().index(i, j);
}

const form_factor::lookup::cross::table_t& form_factor::lookup::cross::raw::get_table() {
    return lookup::detail::get_exv_table_for_settings(
        ff_ax_default_table, lookup::exv::detail::custom_table_ax_raw,
        [](const form_factor::detail::ExvFormFactorSet& set) { return lookup::detail::generate_cross_table<lookup::cross::table_t, RawFormFactorLookup>(set); }
    );
}