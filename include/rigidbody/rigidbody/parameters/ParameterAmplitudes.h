// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/RigidbodyFwd.h>
#include <utility/observer_ptr.h>

namespace ausaxs::rigidbody::parameter {
    /**
     * @brief The maximum magnitude of each parameter component generated per optimisation step.
     *
     * Each amplitude bounds the size of the whole step, not that of its individual coordinates: a translation
     * amplitude of 5 permits a displacement of at most 5 Å, in a uniformly drawn direction.
     */
    struct ParameterAmplitudes {
        double translation = 0;          // Body translation, in Ångström.
        double rotation = 0;             // Body rotation, in radians.
        double symmetry_translation = 0; // Symmetry offset translation, in Ångström.
        double symmetry_rotation = 0;    // Symmetry rotation, in radians.
    };

    /**
     * @brief The body translation amplitude used when none is named.
     */
    constexpr double default_translation() {return 3;}

    /**
     * @brief The body rotation amplitude used when none is named.
     */
    double default_rotation();

    /**
     * @brief The symmetry translation amplitude used when symmetry optimisation is enabled without naming an amplitude.
     */
    double default_symmetry_translation(observer_ptr<const Rigidbody> rigidbody);

    /**
     * @brief The symmetry rotation amplitude used when symmetry optimisation is enabled without naming an amplitude.
     */
    constexpr double default_symmetry_rotation() {return 1.5;}

    /**
     * @brief Every component at its default amplitude, for an optimisation that was not configured by a script.
     */
    ParameterAmplitudes default_amplitudes(observer_ptr<const Rigidbody> rigidbody);
}
