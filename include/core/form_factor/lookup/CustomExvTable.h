// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <form_factor/lookup/FormFactorProduct.h>
#include <form_factor/lookup/ExvFormFactorProduct.h>
#include <form_factor/ExvFormFactor.h>
#include <settings/MoleculeSettings.h>
#include <constants/Constants.h>

namespace ausaxs::form_factor::lookup::exv {
    /**
     * @brief Set custom excluded volume form factor tables for both raw and normalized form factors.
     *        This ensures that raw and normalized form factors stay in sync.
     * 
     * @param set The custom excluded volume form factor set.
     */
    void set_custom_table(const constants::exv::detail::ExvSet& set);

    namespace detail {
        // Shared custom tables for raw form factors
        inline form_factor::lookup::exv::table_t custom_table_xx_raw;
        inline form_factor::lookup::cross::table_t custom_table_ax_raw;

        // Shared custom tables for normalized form factors
        inline form_factor::lookup::exv::table_t custom_table_xx_normalized;
        inline form_factor::lookup::cross::table_t custom_table_ax_normalized;
    }
}