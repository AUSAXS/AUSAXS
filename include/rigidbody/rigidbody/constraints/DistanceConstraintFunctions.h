// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

namespace ausaxs::rigidbody::constraints::functions {
    inline double between_atoms(double offset) {
        return offset*offset*offset*offset*10;
    }

    inline double between_bodies(double offset) {
        return offset*offset*offset*offset*10;
    }

    inline double attractor_repulsor(double offset) {
        return offset*offset*10;
    }
}