// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/DataFwd.h>
#include <data/symmetry/Symmetry.h>
#include <utility/observer_ptr.h>

namespace ausaxs::symmetry {
    enum class type {
        c2,
        c3,
        c4,
        c5,
        c6,
        p2  // dimer with freely-optimizable rotation angle
    };

    /**
     * @brief Apply a given symmetry to a body.
     */
    symmetry::Symmetry get(type t);

    /**
     * @brief Get a predefined symmetry by name.
     */
    type get(std::string_view name);
}