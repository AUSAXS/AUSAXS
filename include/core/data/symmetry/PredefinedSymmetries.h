// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/symmetry/ISymmetry.h>

#include <memory>
#include <string_view>

namespace ausaxs::symmetry {
    enum class type {
        c2,
        c3,
        c4,
        c5,
        c6,
        p2  // dimer with freely-optimizable orientation (PointSymmetry)
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
}
