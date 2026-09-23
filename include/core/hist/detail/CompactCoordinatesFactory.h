// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/DataFwd.h>
#include <form_factor/FormFactorType.h>
#include <form_factor/lookup/FormFactorManager.h>
#include <hist/detail/CompactCoordinates.h>
#include <math/Vector3.h>
#include <utility/Exceptions.h>
#include <utility/observer_ptr.h>

#include <vector>

/**
 * @brief The only construction path for the compact coordinate representations.
 */
namespace ausaxs::hist::detail::factory {
    /**
     * @brief Construct a weight-based representation of @a atoms.
     */
    template<bool variable_bin_width>
    CompactCoordinates<variable_bin_width> construct(const std::vector<data::AtomFF>& atoms) {
        CompactCoordinates<variable_bin_width> c;
        c.fill(atoms);
        return c;
    }

    /**
     * @brief Construct a weight-based representation of every atom in @a molecule.
     */
    template<bool variable_bin_width>
    CompactCoordinates<variable_bin_width> construct_from_atoms(observer_ptr<const data::Molecule> molecule) {
        CompactCoordinates<variable_bin_width> c;
        c.fill_from_atoms(molecule);
        return c;
    }

    /**
     * @brief Construct a weight-based representation of every water in @a molecule.
     */
    template<bool variable_bin_width>
    CompactCoordinates<variable_bin_width> construct_from_waters(observer_ptr<const data::Molecule> molecule) {
        CompactCoordinates<variable_bin_width> c;
        c.fill_from_waters(molecule);
        return c;
    }

    /**
     * @brief Construct a unit-weight representation of @a points.
     *
     * The form factor-aware managers count pairs rather than weigh them, since the form factor amplitudes are applied
     * later, by the intensity calculator. Unit weights make every pair count once.
     */
    template<bool variable_bin_width>
    CompactCoordinates<variable_bin_width> construct_unit_weight(const std::vector<Vector3<double>>& points) {
        CompactCoordinates<variable_bin_width> c;
        c.resize(static_cast<int>(points.size()));
        for (int i = 0; i < c.size(); ++i) {
            c.set_position(i, points[i]);
            c.get_weight(i) = 1;
        }
        return c;
    }

    /**
     * @brief Construct a unit-weight representation of every water in @a molecule. See construct_unit_weight.
     */
    template<bool variable_bin_width>
    CompactCoordinates<variable_bin_width> construct_unit_weight_from_waters(observer_ptr<const data::Molecule> molecule) {
        CompactCoordinates<variable_bin_width> c;
        c.resize(molecule->size_water());
        int i = 0;
        for (const auto& w : molecule->iterate_waters()) {
            c.set_position(i, w.coordinates());
            c.get_weight(i++) = 1;
        }
        return c;
    }

    /**
     * @brief Construct a unit-weight representation of every atom in @a molecule, split by form factor type.
     *        The result is indexed by active form factor index, and types with no atoms get an empty set.
     *        See construct_unit_weight.
     *
     * @throws except::runtime_error if any atom has an UNKNOWN form factor type.
     */
    template<bool variable_bin_width>
    std::vector<CompactCoordinates<variable_bin_width>> construct_unit_weight_by_ff_from_atoms(observer_ptr<const data::Molecule> molecule) {
        auto map = form_factor::manager::get_active_mapping();
        auto active_index = [&map] (const data::AtomFF& a) {
            if (a.form_factor_type() == form_factor::form_factor_t::UNKNOWN) {
                throw except::runtime_error(
                    "factory::construct_unit_weight_by_ff_from_atoms: Attempted to use an atom with UNKNOWN form factor type.\n"
                    "Form factor information is required for the selected excluded volume model."
                );
            }
            return map[static_cast<int>(a.form_factor_type())];
        };

        int n_ff = form_factor::get_active_count();
        std::vector<int> counts(n_ff, 0);
        for (const auto& a : molecule->iterate_atoms()) {++counts[active_index(a)];}

        std::vector<CompactCoordinates<variable_bin_width>> parts(n_ff);
        for (int ff = 0; ff < n_ff; ++ff) {parts[ff].resize(counts[ff]);}

        std::vector<int> filled(n_ff, 0);
        for (const auto& a : molecule->iterate_atoms()) {
            int ff = active_index(a);
            int k = filled[ff]++;
            parts[ff].set_position(k, a.coordinates());
            parts[ff].get_weight(k) = 1;
        }
        return parts;
    }
}
