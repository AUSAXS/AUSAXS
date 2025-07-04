// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/distribution/Distribution2D.h>
#include <hist/distribution/WeightedDistribution2D.h>

namespace ausaxs::hist {
    template <bool UseWeightedContainer>
    struct GenericDistribution2D;

    template <>
    struct GenericDistribution2D<true> {
        using type = WeightedDistribution2D;
    };

    template <>
    struct GenericDistribution2D<false> {
        using type = Distribution2D;
    };
}