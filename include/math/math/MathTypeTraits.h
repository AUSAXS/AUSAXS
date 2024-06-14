#pragma once

#include <type_traits>

namespace detail {
    template<class P>
    struct M {
        operator P const&();
        operator P&&();
    };
    template<class T> constexpr bool is_move_constructible_v = std::is_nothrow_move_constructible_v<T> && !std::is_constructible_v<T, M<T>>;
}
template<class T> constexpr bool supports_nothrow_move_v = detail::is_move_constructible_v<T> && std::is_nothrow_move_assignable_v<T>;