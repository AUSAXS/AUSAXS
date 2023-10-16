#pragma once

#include <cmath>

namespace math {
    namespace fast {
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

    /**
     * @brief Evaluate the exponential function using a Taylor series.
     *        This is probably slower than std::exp, but it is constexpr.
     */
    constexpr double exp(double x) { 
        double sum = 1;
        for (int i = 9; i > 0; --i) {
            sum = 1 + x*sum/i;
        } 
        return sum; 
    } 
}