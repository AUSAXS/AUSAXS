// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/distribution/Distribution3D.h>
#include <hist/distribution/WeightedDistribution3D.h>

namespace ausaxs::hist {
    template <bool UseWeightedContainer, Shape S = Shape::Square>
    struct GenericDistribution3D;

    template <Shape S>
    struct GenericDistribution3D<true, S> {
        using type = WeightedDistribution3D<S>;
    };

    template <Shape S>
    struct GenericDistribution3D<false, S> {
        using type = Distribution3D<S>;
    };
}