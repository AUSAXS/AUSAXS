// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <math/MathFwd.h>

namespace ausaxs::crystal {
    /**
     * @brief A struct to represent a Miller index. 
     */
    struct Miller {
        Miller();
        Miller(int h, int k, int l);

        Vector3<double> normalize() const;
        double distance(const Miller& other) const;
        double length() const;
        double length2() const;
        bool friedel_equivalent(const Miller& other);

        bool operator==(const Miller& other) const;
        bool operator!=(const Miller& other) const;
        Miller operator*(int n) const;

        int h, k, l;
    };
}