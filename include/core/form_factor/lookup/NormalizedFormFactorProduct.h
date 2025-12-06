// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <form_factor/lookup/FormFactorProduct.h>
#include <form_factor/FormFactorType.h>
#include <container/ArrayContainer2D.h>

namespace ausaxs::form_factor {
    using NormalizedFormFactorProduct = FormFactorProduct;

    namespace lookup::atomic {
        using table_t = container::ArrayContainer2D<NormalizedFormFactorProduct, form_factor::get_count(), form_factor::get_count()>;

        namespace normalized {
            /**
             * @brief Get a precalculated atomic form factor product for a given pair of atomic form factors.
             * 
             * @param i The index of the first atomic form factor.
             * @param j The index of the second atomic form factor.
             */
            const NormalizedFormFactorProduct& get_product(unsigned int i, unsigned int j) noexcept;

            /**
             * @brief Get the precalculated atomic form factor product table.
             *        The table is symmetric. 
             */
            const table_t& get_table() noexcept;
        }
    }
}