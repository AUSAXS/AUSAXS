// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <math/Matrix.h>
#include <math/Vector3.h>

namespace ausaxs::transform {
    /**
     * @brief A rigid affine map  v -> rotation*v + translation.
     */
    struct Affine {
        Matrix<double> rotation = Matrix<double>::identity(3);
        Vector3<double> translation{0, 0, 0};

        Vector3<double> operator()(const Vector3<double>& v) const {return rotation*v + translation;}
    };
}
