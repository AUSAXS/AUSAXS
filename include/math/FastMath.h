#pragma once

#include <cmath>

namespace math::fast {
    /**
     * @brief Fast cosine function with a maximum error of 0.001.
     */
    [[maybe_unused]] inline constexpr double cos(double x) noexcept {
        // this is taken directly from https://stackoverflow.com/a/28050328
        constexpr double tp = 1./(2.*M_PI);
        x *= tp;
        x -= 0.25 + std::floor(x + 0.25);
        x *= 16.*(std::abs(x) - 0.5);
        x += .225*x*(std::abs(x) - 1.);
        return x;
    }

    /**
     * @brief Fast sine function with a maximum error of 0.001.
     */
    inline constexpr double sin(double x) noexcept {
        // this is a slightly altered version of https://stackoverflow.com/a/28050328
        constexpr double tp = 1./(2.*M_PI);
        x *= tp;
        x -= 0.5 + std::floor(x);
        x *= 16.*(std::abs(x) - 0.5);
        x += .225*x*(std::abs(x) - 1.);
        return x;
    }
}