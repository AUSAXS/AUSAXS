// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <constants/ConstantsAxes.h>
#include <container/Container2D.h>
#include <hist/distribution/WeightedDistribution2D.h>
#include <utility/TypeTraits.h>

namespace ausaxs::hist {
    /**
     * @brief This is a small wrapper around the Container2D class, indicating that the data
     *        is distributed along the constants::axes::d_vals axis.
     */
    class Distribution2D : public container::Container2D<double> {
        public:
            using Container2D::Container2D;
            explicit Distribution2D(const WeightedDistribution2D& other);
    };
    static_assert(supports_nothrow_move_v<Distribution2D>, "Distribution2D should support nothrow move semantics.");
}