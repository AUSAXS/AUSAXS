// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/RigidbodyFwd.h>
#include <utility/observer_ptr.h>

namespace ausaxs::rigidbody::parameter {
    /**
     * @brief The maximum amplitude of each parameter component generated per optimisation step.
     *
     * A component with a zero amplitude is not generated at all: its field is left unset in the resulting
     * BodyTransformParametersRelative rather than being filled with zeros. There is thus no separate
     * on/off switch; a component is active exactly when its amplitude is non-zero.
     */
    struct ParameterAmplitudes {
        double translation = 0;          // Body translation, in Ångström.
        double rotation = 0;             // Body rotation, in radians.
        double symmetry_translation = 0; // Symmetry offset translation, in Ångström.

        // Reorientation of the symmetry frame. This is in radians for the symmetries parameterised by Euler angles (p2, polyhedral, dihedral),
        // but in raw axis-vector components for the cyclic symmetries, whose rotation angle is instead fixed by their type.
        double symmetry_rotation = 0;
    };

    /**
     * @brief The symmetry translation amplitude used when symmetry optimisation is enabled without naming an amplitude.
     *
     * Scaled by the radius of gyration so that the copies can be displaced across the full extent of the molecule.
     */
    double default_symmetry_translation(observer_ptr<const Rigidbody> rigidbody);

    // The symmetry rotation amplitude used when symmetry optimisation is enabled without naming an amplitude.
    inline constexpr double default_symmetry_rotation = 3;
}
