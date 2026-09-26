// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/DataFwd.h>
#include <hist/detail/CompactCoordinates.h>
#include <math/Vector3.h>
#include <utility/observer_ptr.h>

#include <vector>

/**
 * @brief The only construction path for the compact coordinate representations.
 */
namespace ausaxs::hist::detail::factory {
    /**
     * @brief Construct a weight-based representation of @a atoms.
     */
    inline CompactCoordinates construct(const std::vector<data::AtomFF>& atoms) {
        CompactCoordinates c;
        c.fill(atoms);
        return c;
    }

    /**
     * @brief Construct a weight-based representation of every atom in @a molecule.
     */
    inline CompactCoordinates construct_from_atoms(observer_ptr<const data::Molecule> molecule) {
        CompactCoordinates c;
        c.fill_from_atoms(molecule);
        return c;
    }

    /**
     * @brief Construct a weight-based representation of every water in @a molecule.
     */
    inline CompactCoordinates construct_from_waters(observer_ptr<const data::Molecule> molecule) {
        CompactCoordinates c;
        c.fill_from_waters(molecule);
        return c;
    }

    /**
     * @brief Construct a representation of @a points. A bare point has no weight of its own, so each weighs 1.
     */
    inline CompactCoordinates construct(const std::vector<Vector3<double>>& points) {
        CompactCoordinates c;
        c.resize(static_cast<int>(points.size()));
        for (int i = 0; i < c.size(); ++i) {
            c.set_position(i, points[i]);
            c.get_weight(i) = 1;
        }
        return c;
    }
}
