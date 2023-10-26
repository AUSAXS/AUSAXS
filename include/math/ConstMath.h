#pragma once

#include <math/detail/exp/exp.h>

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
     * @brief Simple constexpr power function.
     *        Only supports integer exponents.  
     */
    constexpr double pow(double x, int n) {
        double result = 1;
        for(int i = 0; i < n; ++i) {
            result *= x;
        }
        if (n < 0) {result = 1./result;}
        return result;
    }

    /**
     * @brief Evaluate the exponential function using IEE4 manipulations, see https://github.com/simonpf/fastexp/tree/master.
     *        This is probably slower than std::exp, but it is constexpr.
     * 
     * @tparam degree The degree of the polynomial fit. This determines the accuracy of the approximation.
     * @param x The argument of the exponential function.
     */
    template<int degree = 5>
    constexpr double exp(double x) {
        return ::detail::fastexp::exp<double, ::detail::fastexp::IEEE, degree>(x);
    }
}