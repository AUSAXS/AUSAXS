// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <form_factor/lookup/CustomExvTable.h>
#include <form_factor/lookup/detail/FormFactorProductBase.h>
#include <form_factor/lookup/detail/LookupHelpers.h>
#include <form_factor/FormFactor.h>
#include <form_factor/NormalizedFormFactor.h>
#include <settings/MoleculeSettings.h>

using namespace ausaxs;
using namespace ausaxs::form_factor;

void form_factor::lookup::exv::set_custom_table(const constants::exv::detail::ExvSet& set) {
    form_factor::detail::ExvFormFactorSet ff(set);

    detail::custom_table_xx_raw = lookup::detail::generate_exv_table<lookup::exv::table_t>(ff);
    detail::custom_table_ax_raw = lookup::detail::generate_cross_table<lookup::cross::table_t, lookup::detail::RawFormFactorLookup>(ff);

    detail::custom_table_xx_normalized = lookup::detail::generate_exv_table<lookup::exv::table_t>(ff);
    detail::custom_table_ax_normalized = lookup::detail::generate_cross_table<lookup::cross::table_t, lookup::detail::NormalizedFormFactorLookup>(ff);

    // ensure that the newly created table will be used in future lookups
    settings::molecule::exv_set = settings::molecule::ExvSet::Custom;
}