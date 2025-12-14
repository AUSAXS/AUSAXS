// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/detail/data/IntrinsicMacros.h>

#include <cstdint>
#include <cstdlib>

namespace ausaxs::hist::detail {
    static inline float squared_dot_product(const float* v1, const float* v2) noexcept {
        float dx = v1[0] - v2[0];
        float dy = v1[1] - v2[1];
        float dz = v1[2] - v2[2];
        return dx*dx + dy*dy + dz*dz;
    }
}

namespace ausaxs::hist::detail {
    // This enum defines the control bits for the SSE4.1 _mm_dp_ps intrinsic
    // the bits are: wzyx | dcba, where the upper nibble selects which components to multiply,
    // and the lower nibble selects where to store the result. This enum therefore defines
    // a multiplication of x, y, z components, and stores the result in the desired location (FIRST, SECOND, ...).
    enum OutputControl : std::int8_t {
        ALL =    0b01111111,
        FIRST =  0b01110001,
        SECOND = 0b01110010,
        THIRD =  0b01110100,
        FOURTH = 0b01111000
    };
}

#if defined __SSE2__
#include <nmmintrin.h>
    namespace ausaxs::hist::detail {
        /**
        * @brief Calculate the squared distance between two CompactCoordinatesXYZW using 128-bit SSE2 instructions.
        */
        static inline __m128 squared_dot_product(const float* v1, const float* v2, OutputControl control) {
            // load data into SSE registers
            __m128 sv1 = _mm_load_ps(v1);
            __m128 sv2 = _mm_load_ps(v2);
            #if defined __SSE4_1__
                // slightly more efficient intrinsic. the weird switch is required with Clang
                __m128 diff = _mm_sub_ps(sv1, sv2);
                switch (control) {
                    case OutputControl::ALL:    return _mm_dp_ps(diff, diff, OutputControl::ALL);
                    case OutputControl::FIRST:  return _mm_dp_ps(diff, diff, OutputControl::FIRST);
                    case OutputControl::SECOND: return _mm_dp_ps(diff, diff, OutputControl::SECOND);
                    case OutputControl::THIRD:  return _mm_dp_ps(diff, diff, OutputControl::THIRD);
                    case OutputControl::FOURTH: return _mm_dp_ps(diff, diff, OutputControl::FOURTH);
                    default: std::abort();
                }
            #else
                sv1[3] = sv2[3] = 0;                        // zero out the weights
                __m128 diff = _mm_sub_ps(sv1, sv2);         // calculate the difference
                __m128 multiplied = _mm_mul_ps(diff, diff); // square the difference
                return _mm_hadd_ps(multiplied, multiplied); // sum the components
            #endif
        }
    }
#endif

#if defined __AVX__
    #include <immintrin.h>
    namespace ausaxs::hist::detail {
        /**
        * @brief Calculate the squared distance between three CompactCoordinatesXYZW using AVX instructions.
        */
        static inline __m256 squared_dot_product(const float* v, const float* v1, const float* v2, OutputControl control) {
            // load data into the 256 bit registers
            __m128 sv = _mm_load_ps(v);
            __m128 sv1 = _mm_load_ps(v1);
            __m128 sv2 = _mm_load_ps(v2);

            __m256 svv = _mm256_broadcast_ps(&sv);      // copy the first 128 bits to the second 128 bits
            __m256 sv12 = _mm256_set_m128(sv2, sv1);    // combine the two 128 bit registers
            __m256 diff = _mm256_sub_ps(svv, sv12);     // calculate the difference
            switch (control) {
                case OutputControl::ALL:    return _mm256_dp_ps(diff, diff, OutputControl::ALL);
                case OutputControl::FIRST:  return _mm256_dp_ps(diff, diff, OutputControl::FIRST);
                case OutputControl::SECOND: return _mm256_dp_ps(diff, diff, OutputControl::SECOND);
                case OutputControl::THIRD:  return _mm256_dp_ps(diff, diff, OutputControl::THIRD);
                case OutputControl::FOURTH: return _mm256_dp_ps(diff, diff, OutputControl::FOURTH);
                default: std::abort();
            }
        }
    }
#endif
