#pragma once

#include <gcem.hpp>

#include <numbers>

namespace ausaxs::gcem {
    namespace fast {
        /**
         * @brief Fast cosine function with a maximum error of 0.001.
         */
        [[maybe_unused]] inline constexpr double cos(double x) noexcept {
            // this is taken directly from https://stackoverflow.com/a/28050328
            constexpr double tp = 1./(2*std::numbers::pi);
            x *= tp;
            x -= 0.25 + ::gcem::floor(x + 0.25);
            x *= 16.*(::gcem::abs(x) - 0.5);
            x += .225*x*(::gcem::abs(x) - 1.);
            return x;
        }

        /**
         * @brief Fast sine function with a maximum error of 0.001.
         */
        inline constexpr double sin(double x) noexcept {
            // this is a slightly altered version of https://stackoverflow.com/a/28050328
            constexpr double tp = 1./(2*std::numbers::pi);
            x *= tp;
            x -= 0.5 + ::gcem::floor(x);
            x *= 16.*(::gcem::abs(x) - 0.5);
            x += .225*x*(::gcem::abs(x) - 1.);
            return x;
        }
    }
}

namespace ausaxs {
    namespace constexpr_math = ::gcem;
}