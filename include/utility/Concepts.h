#pragma once

#include <concepts>

template<typename C>
concept numeric = std::is_arithmetic_v<C>;

template<typename C>
concept indexable = requires(C c, int i) {
    c[i];
};