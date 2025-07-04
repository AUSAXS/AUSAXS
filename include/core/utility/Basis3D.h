// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <math/Vector3.h>

namespace ausaxs {
    class Basis3D {
        public:
            Basis3D() = default;
            Basis3D(const Vector3<double>& x, const Vector3<double>& y, const Vector3<double>& z) : x(x), y(y), z(z) {}
            Vector3<double> x, y, z;
    };
}