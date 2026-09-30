// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <constants/ConstantsAxes.h>
#include <container/ArrayContainer2D.h>
#include <data/DataFwd.h>
#include <form_factor/FormFactorType.h>
#include <form_factor/lookup/FormFactorLookupFwd.h>
#include <utility/observer_ptr.h>

#include <array>
#include <vector>

namespace ausaxs::form_factor::manager {
    namespace detail {
        using profile_t = std::array<double, constants::axes::q_axis.bins>; // A single function evaluated over the default q axis.

        /**
         * @brief The form factor tables of the active form factor set, for the probe selected by settings::scattering::radiation.
         *        All tables are indexed by active slot, not by form_factor_t.
         */
        struct ActiveTables {
            ActiveTables(const std::array<int, form_factor::total_ff_count>& ff_indices, int active_count);
            int active_count;
            std::array<int, form_factor::total_ff_count> ff_indices;
            lookup::table_t raw_exv_table;
            lookup::table_t raw_cross_table;
            lookup::table_t raw_atomic_table;
            lookup::table_t normalized_cross_table;  // Only defined for X-rays, since a neutron form factor may vanish at q = 0.
            lookup::table_t normalized_atomic_table; // Only defined for X-rays, since a neutron form factor may vanish at q = 0.

            /**
             * @brief The amplitude f_i(q) of each slot. The products of these make up raw_atomic_table.
             */
            std::array<profile_t, form_factor::total_ff_count> atomic_profiles{};

            /**
             * @brief The self-term correction s_i(q) - f_i(q)^2 of each slot, where s_i is the scattering of a single group with itself.
             *        raw_atomic_table uses f_i(q)^2 for every pair of the same type, which is only exact for spherically symmetric scatterers. 
             *        This must be added once for every group of the slot, i.e. weighted by the zero-distance bin of its diagonal partial histogram. 
             *        It is only non-zero if self_corrected is true.
             */
            std::array<profile_t, form_factor::total_ff_count> self_correction{};
            bool self_corrected = false;
        };

        /**
         * @brief Activate a custom form factor set.
         *        form_factor_t::OTHER is appended if it is not already present.
         *        With the Fraser-based excluded volume models, form factors without a volume in the current set are removed and treated as OTHER.
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
     *        Requesting the set that is already active is a no-op, so this may be called before every calculation.
     */
    void use_form_factors(const data::Molecule& molecule);

    /**
     * @brief Rebuild the active product tables in-place, preserving the current form factor selection.
     *        Called whenever the EXV method or parameter set changes, since this may also change which form factors are available (see detail::use_form_factors).
     */
    void rebuild();
}
