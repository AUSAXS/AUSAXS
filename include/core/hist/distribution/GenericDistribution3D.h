// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/distribution/Distribution3D.h>
#include <hist/distribution/WeightedDistribution3D.h>

namespace ausaxs::hist {
    template <bool UseWeightedContainer>
    struct GenericDistribution3D;

    template <>
    struct GenericDistribution3D<true> {
        using type = WeightedDistribution3D;
    };

    template <>
    struct GenericDistribution3D<false> {
        using type = Distribution3D;
    };
}