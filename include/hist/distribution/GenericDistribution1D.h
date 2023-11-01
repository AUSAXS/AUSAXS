#pragma once

#include <hist/distribution/Distribution1D.h>
#include <hist/distribution/WeightedDistribution1D.h>

namespace hist {
    template <bool UseWeightedContainer, typename T>
    struct GenericDistribution1D;

    template <typename T>
    struct GenericDistribution1D<true, T> {
        using type = WeightedDistribution1D<T>;
    };

    template <typename T>
    struct GenericDistribution1D<false, T> {
        using type = Distribution1D<T>;
    };
}