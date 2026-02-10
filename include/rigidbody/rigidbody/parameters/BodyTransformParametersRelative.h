// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <math/Vector3.h>
#include <data/symmetry/Symmetry.h>

#include <vector>
#include <optional>

namespace ausaxs::rigidbody::parameter {
    /**
     * @brief A relative (delta) transformation to be applied to a body or rigid group.
     * 
     * This represents a change in position/orientation, not an absolute state.
     * The rotation is an axis-angle representation of the delta rotation.
     * The translation is the delta translation vector.
     */
    struct BodyTransformParametersRelative {
        BodyTransformParametersRelative() = default;
        BodyTransformParametersRelative(const Vector3<double>& translation, const Vector3<double>& rotation, std::vector<symmetry::Symmetry>&& symmetry_pars = {})
            : translation(translation), rotation(rotation), symmetry_pars(std::move(symmetry_pars)) 
        {}
        std::optional<Vector3<double>> translation;
        std::optional<Vector3<double>> rotation;
        std::optional<std::vector<symmetry::Symmetry>> symmetry_pars;

        bool operator==(const BodyTransformParametersRelative&) const = default;
    };
}
