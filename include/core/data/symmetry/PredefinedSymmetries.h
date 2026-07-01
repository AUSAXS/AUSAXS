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
        icosahedral,// chiral icosahedral rotation group (PolyhedralSymmetry)
        d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12 // chiral dihedral rotation groups (DihedralSymmetry)
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
     * A plain name (e.g. "c3", "p2", "t", "d3") maps to the corresponding predefined symmetry.
     * A hyphenated name builds a nested CompositeSymmetry: the first part is the inner
     * symmetry, the remainder (recursively parsed) the outer one. For example "p2-c3" nests
     * a p2 dimer inside an outer c3, and "c2-c2-c3" nests left-to-right as c2-(c2-c3).
     *
     * The one exception is a bare "c2-cN" / "cN-c2" pair (either order, N>=2): a C2 and a single
     * CN sharing a centre with perpendicular axes is the dihedral group D_N, so it is built as a
     * DihedralSymmetry rather than a generic composite. This gives the correct 2N-copy point group
     * (instead of the over-parameterised composite) and lets the pair-distance schedule exploit the
     * full group symmetry. "c2-c2" is D2 accordingly.
     */
    std::unique_ptr<ISymmetry> create(std::string_view name);
}
