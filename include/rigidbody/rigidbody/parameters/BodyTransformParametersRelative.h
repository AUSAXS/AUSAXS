// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <math/Vector3.h>
#include <data/symmetry/ISymmetry.h>

#include <vector>
#include <optional>

namespace ausaxs::rigidbody::parameter {
    /**
     * @brief A relative (delta) transformation to be applied to a body or rigid group.
     */
    struct BodyTransformParametersRelative {
        BodyTransformParametersRelative();
        BodyTransformParametersRelative(
            const Vector3<double>& translation, 
            const Vector3<double>& rotation, 
            std::vector<std::unique_ptr<symmetry::ISymmetry>> symmetry_pars = {}
        );
        
        std::optional<Vector3<double>> translation;
        std::optional<Vector3<double>> rotation;
        std::optional<std::vector<std::unique_ptr<symmetry::ISymmetry>>> symmetry_pars;
        bool operator==(const BodyTransformParametersRelative&) const;
    };
}
