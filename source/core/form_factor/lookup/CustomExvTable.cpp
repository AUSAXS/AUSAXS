// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include "settings/MoleculeSettings.h"
#include <form_factor/lookup/CustomExvTable.h>
#include <form_factor/lookup/detail/FormFactorProductBase.h>
#include <form_factor/FormFactor.h>
#include <form_factor/NormalizedFormFactor.h>
#include <settings/ExvSettings.h>

using namespace ausaxs;
using namespace ausaxs::form_factor;

namespace {
    struct RawFormFactorLookup {
        static constexpr const FormFactor& get(form_factor_t type) {
            return lookup::atomic::raw::get(type);
        }
    };

    struct NormalizedFormFactorLookup {
        static constexpr const NormalizedFormFactor& get(form_factor_t type) {
            return lookup::atomic::normalized::get(type);
        }
    };
}

void form_factor::lookup::exv::set_custom_table(const constants::exv::detail::ExvSet& set) {
    form_factor::detail::ExvFormFactorSet ff(set);

    detail::custom_table_xx_raw = lookup::detail::generate_exv_table<lookup::exv::table_t>(ff);
    detail::custom_table_ax_raw = lookup::detail::generate_cross_table<lookup::cross::table_t, RawFormFactorLookup>(ff);

    detail::custom_table_xx_normalized = lookup::detail::generate_exv_table<lookup::exv::table_t>(ff);
    detail::custom_table_ax_normalized = lookup::detail::generate_cross_table<lookup::cross::table_t, NormalizedFormFactorLookup>(ff);

    // ensure that the newly created table will be used in future lookups
    settings::molecule::exv_set = settings::molecule::ExvSet::Custom;
}