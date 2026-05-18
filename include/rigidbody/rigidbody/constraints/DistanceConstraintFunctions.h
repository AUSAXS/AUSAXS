// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

// Penalty functions for the various distance constraints. Each maps the deviation of a constraint
// from its target distance (`offset`) to a non-negative penalty added to the optimisation score.
namespace ausaxs::rigidbody::constraints::functions {
    /**
     * @brief Penalty for a distance constraint between two individual atoms. Quartic in the offset.
     */
    inline double between_atoms(double offset) {
        return offset*offset*offset*offset*10;
    }

    /**
     * @brief Penalty for a distance constraint between two bodies. Quartic in the offset.
     */
    inline double between_bodies(double offset) {
        return offset*offset*offset*offset*10;
    }

    /**
     * @brief Penalty applied by attractor and repulsor constraints. Quadratic in the offset.
     */
    inline double attractor_repulsor(double offset) {
        return offset*offset*10;
    }
}