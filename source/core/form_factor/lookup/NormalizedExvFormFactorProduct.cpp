// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <form_factor/lookup/NormalizedExvFormFactorProduct.h>
#include <form_factor/lookup/CustomExvTable.h>
#include <form_factor/lookup/detail/FormFactorProductBase.h>
#include <form_factor/lookup/detail/LookupHelpers.h>
#include <form_factor/NormalizedFormFactor.h>
#include <settings/MoleculeSettings.h>
#include <constants/Constants.h>

using namespace ausaxs;
using namespace ausaxs::form_factor;

namespace {
    auto ff_xx_default_table = lookup::detail::generate_exv_table<lookup::exv::table_t>(lookup::exv::standard);
    auto ff_ax_default_table = lookup::detail::generate_cross_table<lookup::cross::table_t, lookup::detail::NormalizedFormFactorLookup>(lookup::exv::standard);
}

const NormalizedFormFactorProduct& form_factor::lookup::exv::normalized::get_product(unsigned int i, unsigned int j) {
    return get_table().index(i, j);
}

const form_factor::lookup::exv::table_t& form_factor::lookup::exv::normalized::get_table() {
    return lookup::detail::get_exv_table_for_settings(
        ff_xx_default_table, lookup::exv::detail::custom_table_xx_normalized,
        [](const form_factor::detail::ExvFormFactorSet& set) { return lookup::detail::generate_exv_table<lookup::exv::table_t>(set); }
    );
}

const NormalizedFormFactorProduct& form_factor::lookup::cross::normalized::get_product(unsigned int i, unsigned int j) {
    return get_table().index(i, j);
}

const form_factor::lookup::cross::table_t& form_factor::lookup::cross::normalized::get_table() {
    return lookup::detail::get_exv_table_for_settings(
        ff_ax_default_table, lookup::exv::detail::custom_table_ax_normalized,
        [](const form_factor::detail::ExvFormFactorSet& set) { return lookup::detail::generate_cross_table<lookup::cross::table_t, lookup::detail::NormalizedFormFactorLookup>(set); }
    );
}