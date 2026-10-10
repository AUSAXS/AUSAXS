// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <constants/ConstantsAxes.h>
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
         *        The product tables hold the cross-correlation form factors of two distinct scatterers.
         *        A scatterer correlated with itself instead uses the self-correlation form factor of its slot.
         */
        struct ActiveTables {
            ActiveTables(const std::array<int, form_factor::total_ff_count>& ff_indices, int active_count);
            int active_count;
            std::array<int, form_factor::total_ff_count> ff_indices;
            lookup::table_t raw_exv_table;
            lookup::table_t raw_cross_table;
            lookup::table_t raw_atomic_table;
            std::vector<profile_t> raw_self_table; // One entry per active slot.
        };

        /**
         * @brief Activate a custom form factor set.
         *        form_factor_t::OTHER is appended if it is not already present.
         *        With the Fraser-based excluded volume models, form factors without a volume in the current set are removed and treated as OTHER.
         *        Throws if the selection, including OTHER, exceeds settings::form_factor::max_types.
         */
        void use_form_factors(std::vector<int> ff_indices);
    }

    /**
     * @brief Evaluate the amplitude f(q) of a single form factor over the default q axis, for the probe selected by settings::scattering::radiation.
     */
    detail::profile_t evaluate_amplitude(form_factor_t type);

    /**
     * @brief Evaluate the self-correlation form factor of a single form factor over the default q axis, for the probe selected by settings::scattering::radiation.
     *        This is the form factor of a scatterer correlated with itself, which is only the squared amplitude for spherically symmetric scatterers.
     */
    detail::profile_t evaluate_self(form_factor_t type);

    /**
     * @brief Get the currently active form factor product tables. 
     *        Throws if no form factor selection has been made yet; see use_form_factors.
     */
    observer_ptr<const detail::ActiveTables> get_active_product_tables();

    /**
     * @brief Get a mapping from form_factor_t enum index to active slot index.
     *        All form factors not in the active set are mapped to OTHER. 
     */
    std::vector<int> get_active_mapping();

    /**
     * @brief Determine the most appropriate form factor set for the given molecule and activate it. 
     *        The most abundant types get a slot of their own, up to settings::form_factor::max_types slots in total. Types rarer than
     *        settings::form_factor::min_fraction of the atoms are folded onto OTHER regardless.
     *        Requesting the set that is already active is a no-op, so this may be called before every calculation.
     */
    void use_form_factors(const data::Molecule& molecule);

    /**
     * @brief Rebuild the active product tables in-place, preserving the current form factor selection.
     *        Called whenever the EXV method or parameter set changes, since this may also change which form factors are available (see detail::use_form_factors).
     */
    void rebuild();
}
