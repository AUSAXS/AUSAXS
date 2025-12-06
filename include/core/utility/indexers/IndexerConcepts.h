#pragma once

#include <concepts>

namespace ausaxs::utility::indexer {
    template<class T>
    concept ValidIndexable1D = requires(T a) {
        { a.data } ;
        { a.N } -> std::convertible_to<int>;
    };

    template<class T>
    concept ValidIndexable2D = requires(T a) {
        { a.M } -> std::convertible_to<int>;
    } && ValidIndexable1D<T>;

    template<class T>
    concept ValidIndexable3D = requires(T a) {
        { a.L } -> std::convertible_to<int>;
    } && ValidIndexable2D<T>;
}