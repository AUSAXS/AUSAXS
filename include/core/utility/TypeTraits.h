#pragma once

#include <math/MathTypeTraits.h>

#include <string_view>

namespace detail {
    template<typename T>
    concept string_type = std::same_as<T, std::string> || std::same_as<T, std::string_view> || std::same_as<T, const char*>;

}