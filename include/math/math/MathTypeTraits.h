#pragma once

#include <type_traits>

namespace ausaxs::detail {
    #if defined(_MSC_VER)
        template<class T> constexpr bool is_move_constructible_v = std::is_nothrow_move_constructible_v<T>;
    #else
        // taken from https://stackoverflow.com/a/51912859
        template<class P>
        struct M {
            operator P const&();
            operator P&&();
        };
        template<class T> constexpr bool is_move_constructible_v = std::is_nothrow_move_constructible_v<T> && !std::is_constructible_v<T, M<T>>;
    #endif
}
template<class T> constexpr bool supports_nothrow_move_v = ausaxs::detail::is_move_constructible_v<T> && std::is_nothrow_move_assignable_v<T>;