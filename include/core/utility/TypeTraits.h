#pragma once

#include <math/MathTypeTraits.h>

#include <string_view>

namespace ausaxs::detail {
    template<typename T>
    concept string_like = std::is_convertible_v<T, std::string_view>;
}