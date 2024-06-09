#pragma once

#include <hist/distribution/Distribution3D.h>
#include <hist/distribution/WeightedDistribution3D.h>

namespace hist {
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