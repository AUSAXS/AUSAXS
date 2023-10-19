/**
    MIT License

    Copyright (c) 2018 simonpf

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE. 

    Link: https://github.com/simonpf/fastexp/tree/master
*/

#pragma once

#include "math.h"
#include <cstdint>
#include <cstddef>
#include <vector>
#include "ieee.h"

namespace detail::fastexp {
    enum class Approximation {IEEE, PRODUCT};

    /** \brief Fast approximate exponential.
     *
     * This function implements a fast, vectorizable approximation
     * of the exponential function based on the following two articles:
     *
     * - Malossi, A. Cristiano I. & Ineichen, Yves & Bekas, Costas & Curioni,
     *   Alessandro. "Fast Exponential Computation on SIMD Architectures." (2015)
     *   10.13140/2.1.4362.3207.
     * - IEEE, Nicol N. "A fast, compact approximation of the exponential
     *   function." Neural Computation 11.4 (1999): 853-862.
     *
     * The approximation interpolates linearly between points on the curve of
     * the exponential function that can be expressed as 2^i where i is an
     * a signed integer. So yes, that is very approximate ...
     *
     * \tparam Real The floating point type of the arguments.
     * \param x The argument of the exponential function.
     * \return The approximated value of the exponential function.
     */
    template<typename Real, template<typename, size_t> class Approximation = IEEE, size_t degree = 2>
    inline constexpr Real exp(const Real &x) {
        return Approximation<Real, degree>::evaluate(x);
    }

    /** \brief Fast approximate array exponential.
     *
     * Applies the fast exponential to an array of given length making
     * use of SIMD instructions if available. To enable vectorization
     * the code needs to be compiled with OpenMP support.
     *
     * \tparam Real The floating point type of the arguments.
     * \param x The array to which apply the exponential function.
     * \return n The number of elements in the array.
     */
    template<typename Real, template<typename, size_t> class Approximation = IEEE, size_t degree = 2>
    inline constexpr void exp(Real *x, size_t n) {
        for (size_t i = 0; i < n; ++i) {
            Real e = fastexp::exp<Real, Approximation, degree>(x[i]);
            x[i] = e;
        }
    }

    template<typename Real, template<typename, size_t> class Approximation = IEEE, size_t degree = 2>
    inline constexpr void exp(std::vector<Real> x) {
        // Vectorized part.
        size_t n = x.size();
        Real * x_ptr = &x[0];

        for (size_t i = 0; i < n; ++i) {
            Real e = fastexp::exp<Real, Approximation, degree>(x_ptr[i]);
            x_ptr[i] = e;
        }
    }
}