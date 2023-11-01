#pragma once

#include <hist/distribution/Distribution2D.h>
#include <hist/distribution/WeightedDistribution2D.h>

namespace hist {
    template <bool UseWeightedContainer, typename T>
    struct GenericDistribution2D;

    template <typename T>
    struct GenericDistribution2D<true, T> {
        using type = WeightedDistribution2D<T>;
    };

    template <typename T>
    struct GenericDistribution2D<false, T> {
        using type = Distribution2D<T>;
    };
}