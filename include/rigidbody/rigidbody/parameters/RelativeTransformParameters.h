// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <math/Vector3.h>
#include <data/symmetry/Symmetry.h>

#include <vector>

namespace ausaxs::rigidbody::parameter {
    /**
     * @brief A relative (delta) transformation to be applied to a body or rigid group.
     * 
     * This represents a change in position/orientation, not an absolute state.
     * The rotation is an axis-angle representation of the delta rotation.
     * The translation is the delta translation vector.
     */
    struct RelativeTransformParameters {
        Vector3<double> translation = {0, 0, 0};
        Vector3<double> rotation = {0, 0, 0};
        std::vector<symmetry::Symmetry> symmetry_pars;

        bool operator==(const RelativeTransformParameters&) const = default;
    };
}
