// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/DataFwd.h>
#include <form_factor/FormFactorType.h>
#include <form_factor/lookup/FormFactorManager.h>
#include <hist/detail/CompactCoordinatesFactory.h>
#include <utility/Exceptions.h>
#include <utility/observer_ptr.h>

#include <vector>

/**
 * @brief The construction of the compact coordinate representations for the form factor-resolved histograms.
 */
namespace ausaxs::hist::detail::factory {
    /**
     * @brief Construct a weight-based representation of every atom in @a molecule, split by form factor type.
     *        The result is indexed by active form factor index, and types with no atoms get an empty set.
     */
    inline std::vector<CompactCoordinates> construct_by_ff_from_atoms(observer_ptr<const data::Molecule> molecule) {
        auto map = form_factor::manager::get_active_mapping();
        auto active_index = [&map] (const data::AtomFF& a) {
            if (a.form_factor_type() == form_factor::form_factor_t::UNKNOWN) {
                throw except::runtime_error(
                    "factory::construct_by_ff_from_atoms: Attempted to use an atom with UNKNOWN form factor type.\n"
                    "Form factor information is required for the selected excluded volume model."
                );
            }
            return map[static_cast<int>(a.form_factor_type())];
        };

        int n_ff = form_factor::get_active_count();
        std::vector<int> counts(n_ff, 0);
        for (const auto& a : molecule->iterate_atoms()) {++counts[active_index(a)];}

        std::vector<CompactCoordinates> parts(n_ff);
        for (int ff = 0; ff < n_ff; ++ff) {parts[ff].resize(counts[ff]);}

        std::vector<int> filled(n_ff, 0);
        for (const auto& a : molecule->iterate_atoms()) {
            int ff = active_index(a);
            int k = filled[ff]++;
            parts[ff].set_position(k, a.coordinates());
            parts[ff].get_weight(k) = static_cast<float>(a.weight());
        }
        return parts;
    }
}
