#pragma once

#include <type_traits>

namespace crystal {
    /**
     * @brief A struct to represent a Miller index. 
     * 
     * @tparam T The type of the Miller index. In most cases this should be left as the default. 
     */
    template<typename T = unsigned, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
    struct Miller {
        Miller() = default;
        Miller(T h, T k, T l) : h(h), k(k), l(l) {}
        T h, k, l;
    };
}