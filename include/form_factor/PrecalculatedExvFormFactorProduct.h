#pragma once

#include <container/ContainerFwd.h>
#include <constants/Constants.h>
#include <form_factor/PrecalculatedFormFactorProduct.h>
#include <form_factor/ExvFormFactor.h>
#include <container/ArrayContainer2D.h>

#include <vector>

namespace form_factor::storage::exv {
    /**
     * @brief Get a precalculated excluded-volume form factor product for a given pair of excluded volume form factors.
     * 
     * @param i The index of the first excluded volume form factor.
     * @param j The index of the second excluded volume form factor.
     */
    const PrecalculatedFormFactorProduct& get_precalculated_form_factor_product(unsigned int i, unsigned int j) noexcept;

    /**
     * @brief Get the precalculated excluded-volume form factor product table.
     *        The table is symmetric. 
     */
    const container::ArrayContainer2D<PrecalculatedFormFactorProduct, form_factor::get_count_without_excluded_volume(), form_factor::get_count_without_excluded_volume()>& get_precalculated_form_factor_table() noexcept;
}

namespace form_factor::storage::cross {
    /**
     * @brief Get a precalculated form factor product for a given pair of excluded volume and atomic form factors. 
     * 
     * @param i The index of the *atomic* form factor.
     * @param j The index of the *excluded volume* form factor. 
     */
    const PrecalculatedFormFactorProduct& get_precalculated_form_factor_product(unsigned int i, unsigned int j) noexcept;

    /**
     * @brief Get the precalculated form factor product table for a given pair of excluded volume and atomic form factors.
     *        The first index is the atomic form factor, and the second index is the excluded volume form factor.
     */
    const container::ArrayContainer2D<PrecalculatedFormFactorProduct, form_factor::get_count_without_excluded_volume(), form_factor::get_count_without_excluded_volume()>& get_precalculated_form_factor_table() noexcept;
}