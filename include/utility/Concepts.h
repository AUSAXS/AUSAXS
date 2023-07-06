#pragma once

#include <concepts>

template<typename C>
concept numeric = std::is_arithmetic_v<C>;

// container_type requires both operator[] and size()
template<typename C>
concept container_type = requires(C c, int i) {
    {c[i]};
    {c.size()};
};

template<typename C>
concept indexable = requires(C c, int i) {
    c[i];
};