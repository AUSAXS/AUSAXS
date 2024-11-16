#pragma once

#include <hist/distribution/Distribution1D.h>
#include <hist/distribution/WeightedDistribution1D.h>

namespace ausaxs::hist {
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