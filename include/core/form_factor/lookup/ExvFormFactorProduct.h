// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <container/ContainerFwd.h>
#include <constants/Constants.h>
#include <form_factor/lookup/FormFactorProduct.h>
#include <form_factor/ExvTable.h>
#include <container/ArrayContainer2D.h>

namespace ausaxs::form_factor::lookup::exv {
    using table_t = container::ArrayContainer2D<FormFactorProduct, form_factor::get_count_without_excluded_volume(), form_factor::get_count_without_excluded_volume()>;

    namespace raw {
        /**
        * @brief Get a precalculated excluded-volume form factor product for a given pair of excluded volume form factors.
        * 
        * @param i The index of the first excluded volume form factor.
        * @param j The index of the second excluded volume form factor.
        */
        const FormFactorProduct& get_product(unsigned int i, unsigned int j);

        /**
        * @brief Get the precalculated excluded-volume form factor product table.
        *        The table is symmetric. 
        */
        const table_t& get_table();
    }
}

namespace ausaxs::form_factor::lookup::cross {
    using table_t = container::ArrayContainer2D<FormFactorProduct, form_factor::get_count_without_excluded_volume(), form_factor::get_count_without_excluded_volume()>;

    namespace raw {
        /**
        * @brief Get a precalculated form factor product for a given pair of excluded volume and atomic form factors. 
        * 
        * @param i The index of the *atomic* form factor.
        * @param j The index of the *excluded volume* form factor. 
        */
        const FormFactorProduct& get_product(unsigned int i, unsigned int j);

        /**
        * @brief Get the precalculated form factor product table for a given pair of excluded volume and atomic form factors.
        *        The first index is the atomic form factor, and the second index is the excluded volume form factor.
        */
        const table_t& get_table();
    }
}

namespace ausaxs::form_factor::lookup::exv::raw::detail {
    void set_custom_table(const constants::exv::detail::ExvSet& set);
}