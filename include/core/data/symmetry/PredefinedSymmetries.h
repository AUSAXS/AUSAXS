// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/symmetry/ISymmetry.h>

#include <memory>
#include <string_view>

namespace ausaxs::symmetry {
    enum class type {
        c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, // cyclic symmetries
        p2, // dimer with freely-optimizable orientation (PointSymmetry)
        tetrahedral,// chiral tetrahedral rotation group (PolyhedralSymmetry)
        octahedral, // chiral octahedral rotation group (PolyhedralSymmetry)
        icosahedral // chiral icosahedral rotation group (PolyhedralSymmetry)
    };

    /**
     * @brief Create a predefined symmetry object of the given type.
     *        Returns a Symmetry for c-type; returns a PointSymmetry for p-type.
     */
    std::unique_ptr<ISymmetry> get(type t);

    /**
     * @brief Get a predefined symmetry enum value by name string.
     */
    type get(std::string_view name);

    /**
     * @brief Build a symmetry object from a name string.
     *
     * A plain name (e.g. "c3", "p2", "t") maps to the corresponding predefined symmetry.
     * A hyphenated name builds a nested CompositeSymmetry: the first part is the inner
     * symmetry, the remainder (recursively parsed) the outer one. For example "p2-c3" nests
     * a p2 dimer inside an outer c3, and "c2-c2-c3" nests left-to-right as c2-(c2-c3).
     */
    std::unique_ptr<ISymmetry> create(std::string_view name);
}
