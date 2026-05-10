// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/DataFwd.h>
#include <container/ArrayContainer2D.h>
#include <form_factor/lookup/FormFactorLookupFwd.h>
#include <utility/observer_ptr.h>

#include <array>
#include <vector>

namespace ausaxs::form_factor::manager {
    namespace detail {
        struct ActiveTables {
            ActiveTables(std::array<int, settings::form_factor::max_ff_types>&& ff_indices, unsigned int active_count);
            unsigned int active_count;
            std::array<int, settings::form_factor::max_ff_types> ff_indices;
            lookup::table_t    raw_exv_table;
            lookup::table_t  raw_cross_table;
            lookup::table_t raw_atomic_table;
            lookup::table_t  normalized_cross_table;
            lookup::table_t normalized_atomic_table;
        };

        /**
         * @brief Activate a custom form factor set.
         */
        void use_form_factors(std::vector<int> ff_indices);
    }

    /**
     * @brief Get the currently active form factor product tables. 
     */
    observer_ptr<const detail::ActiveTables> get_active_product_tables() noexcept;

    /**
     * @brief Get a mapping from form_factor_t enum index to active slot index.
     *        All form factors not in the active set are mapped to OTHER. 
     */
    std::vector<int> get_active_mapping();

    /**
     * @brief Determine the most appropriate form factor set for the given molecule and activate it. 
     */
    void use_form_factors(data::Molecule& molecule);

    /**
     * @brief Rebuild the active product tables in-place, preserving the current form factor selection.
     *        Should be called whenever the EXV parameter set changes.
     */
    void rebuild();
}

namespace ausaxs::form_factor {
    /**
     * @brief Get the number of active form factors. 
     *        This is often smaller than settings::form_factor::max_ff_types when using custom form factors. 
     */
    unsigned int get_active_count() noexcept;
}