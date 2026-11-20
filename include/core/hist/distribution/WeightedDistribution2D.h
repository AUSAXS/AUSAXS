// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <constants/ConstantsAxes.h>
#include <container/Container2D.h>
#include <hist/distribution/detail/WeightedEntry.h>
#include <settings/InternalState.h>
#include <utility/TypeTraits.h>

#include <cmath>

namespace ausaxs::hist {
    /**
     * @brief This is a small wrapper around the Container2D class. Anything added to this
     *        distribution will be tracked by the WeightedDistribution class, which may add
     *        a significant overhead compared to a pure Distribution2D class.
     */
    class WeightedDistribution2D : public container::Container2D<detail::WeightedEntry> {
        public:
            using Container2D::Container2D;
    };
    static_assert(supports_nothrow_move_v<WeightedDistribution2D>, "WeightedDistribution2D should support nothrow move semantics.");
}