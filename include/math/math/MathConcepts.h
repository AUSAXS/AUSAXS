#pragma once

#include <concepts>

template<typename C>
concept numeric = std::is_arithmetic_v<C>;