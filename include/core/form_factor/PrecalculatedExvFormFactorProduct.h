// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <container/ContainerFwd.h>
#include <constants/Constants.h>
#include <form_factor/PrecalculatedFormFactorProduct.h>
#include <form_factor/ExvFormFactor.h>
#include <container/ArrayContainer2D.h>

namespace ausaxs::form_factor::storage::exv {
    using table_t = container::ArrayContainer2D<PrecalculatedFormFactorProduct, form_factor::get_count_without_excluded_volume(), form_factor::get_count_without_excluded_volume()>;
    /**
     * @brief Get a precalculated excluded-volume form factor product for a given pair of excluded volume form factors.
     * 
     * @param i The index of the first excluded volume form factor.
     * @param j The index of the second excluded volume form factor.
     */
    const PrecalculatedFormFactorProduct& get_precalculated_form_factor_product(unsigned int i, unsigned int j);

    /**
     * @brief Get the precalculated excluded-volume form factor product table.
     *        The table is symmetric. 
     */
    const table_t& get_precalculated_form_factor_table();
}

namespace ausaxs::form_factor::storage::cross {
    using table_t = container::ArrayContainer2D<PrecalculatedFormFactorProduct, form_factor::get_count_without_excluded_volume(), form_factor::get_count_without_excluded_volume()>;

    /**
     * @brief Get a precalculated form factor product for a given pair of excluded volume and atomic form factors. 
     * 
     * @param i The index of the *atomic* form factor.
     * @param j The index of the *excluded volume* form factor. 
     */
    const PrecalculatedFormFactorProduct& get_precalculated_form_factor_product(unsigned int i, unsigned int j);

    /**
     * @brief Get the precalculated form factor product table for a given pair of excluded volume and atomic form factors.
     *        The first index is the atomic form factor, and the second index is the excluded volume form factor.
     */
    const table_t& get_precalculated_form_factor_table();
}

namespace ausaxs::form_factor::storage::detail {
    void set_custom_displaced_volume_table(const constants::displaced_volume::detail::DisplacedVolumeSet& set);
}