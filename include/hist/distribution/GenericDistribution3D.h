#pragma once

#include <hist/distribution/Distribution2D.h>
#include <hist/distribution/WeightedDistribution2D.h>

namespace hist {
    template <bool UseWeightedContainer, typename T>
    struct GenericDistribution3D;

    template <typename T>
    struct GenericDistribution3D<true, T> {
        using type = WeightedDistribution3D<T>;
    };

    template <typename T>
    struct GenericDistribution3D<false, T> {
        using type = Distribution3D<T>;
    };
}