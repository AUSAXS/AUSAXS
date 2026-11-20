// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <constants/ConstantsAxes.h>
#include <container/Container3D.h>
#include <hist/distribution/DistributionFwd.h>
#include <hist/distribution/detail/WeightedEntry.h>
#include <settings/InternalState.h>
#include <utility/TypeTraits.h>

namespace ausaxs::hist {
    /**
     * @brief This is a small wrapper around the Container3D class. Anything added to this
     *        distribution will be tracked by the WeightedDistribution class, which may add
     *        a significant overhead compared to a pure Distribution3D class.
     *
     * @tparam S Shape::Triangular stores one histogram per unordered pair (x, y) only; see Distribution3D.
     */
    template<Shape S>
    class WeightedDistribution3D : public container::Container3D<detail::WeightedEntry, S> {
        public:
            WeightedDistribution3D() = default;
            WeightedDistribution3D(int width, int height, int depth) : container::Container3D<detail::WeightedEntry, S>(width, height, depth) {}
            WeightedDistribution3D(int width, int height, int depth, const detail::WeightedEntry& value) : container::Container3D<detail::WeightedEntry, S>(width, height, depth, value) {}

            using container::Container3D<detail::WeightedEntry, S>::index;
            using container::Container3D<detail::WeightedEntry, S>::linear_index;
    };
    static_assert(supports_nothrow_move_v<WeightedDistribution3D<Shape::Square>>, "WeightedDistribution3D should support nothrow move semantics.");
    static_assert(supports_nothrow_move_v<WeightedDistribution3D<Shape::Triangular>>, "WeightedDistribution3D should support nothrow move semantics.");
}