// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/DataFwd.h>
#include <hist/detail/CompactCoordinates.h>
#include <hist/detail/CompactCoordinatesFF.h>
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
     * @brief Construct a form-factor-based representation of @a atoms.
     */
    template<bool variable_bin_width>
    CompactCoordinatesFF<variable_bin_width> construct_ff(const std::vector<data::AtomFF>& atoms) {
        CompactCoordinatesFF<variable_bin_width> c;
        c.fill(atoms);
        c.setup();
        return c;
    }

    /**
     * @brief Construct a form-factor-based representation of every atom in @a molecule.
     */
    template<bool variable_bin_width>
    CompactCoordinatesFF<variable_bin_width> construct_ff_from_atoms(observer_ptr<const data::Molecule> molecule) {
        CompactCoordinatesFF<variable_bin_width> c;
        c.fill_from_atoms(molecule);
        c.setup();
        return c;
    }

    /**
     * @brief Construct a form-factor-based representation of every water in @a molecule.
     */
    template<bool variable_bin_width>
    CompactCoordinatesFF<variable_bin_width> construct_ff_from_waters(observer_ptr<const data::Molecule> molecule) {
        CompactCoordinatesFF<variable_bin_width> c;
        c.fill_from_waters(molecule);
        c.setup();
        return c;
    }
}
