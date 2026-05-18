// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/distribution/Distribution1D.h>
#include <hist/distribution/WeightedDistribution1D.h>

namespace ausaxs::hist {
    /**
     * @brief Trait selecting the 1D distance distribution type to use.
     *
     * Resolves to WeightedDistribution1D when @p UseWeightedContainer is true, and to the plain
     * Distribution1D otherwise. This lets histogram code be templated on whether the per-bin mean
     * distance is tracked (weighted bins) or not (unweighted bins).
     */
    template <bool UseWeightedContainer>
    struct GenericDistribution1D;

    template<>
    struct GenericDistribution1D<true> {
        using type = WeightedDistribution1D;
    };

    template<>
    struct GenericDistribution1D<false> {
        using type = Distribution1D;
    };
}