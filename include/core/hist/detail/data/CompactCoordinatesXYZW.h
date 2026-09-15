// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

/**
 * @brief This file contains custom vector instructions for efficient scattering calculations.
 *        The implementation is specialized for generic systems defined by separate x, y, z and w arrays.
 */

#pragma once

#include <constants/Constants.h>
#include <hist/detail/data/IntrinsicHelpers.h>
#include <hist/detail/data/IntrinsicMacros.h>
#include <hist/detail/data/WidthControllers.h>
#include <settings/InternalState.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <type_traits>

namespace ausaxs::hist::detail::xyzw {
    struct EvaluatedResult {
        float distance;         // The raw distance
        int32_t distance_bin;   // The distance bin index
        float weight;           // The combined weight
    };

    // same as above, except it does not provide the exact distance
    struct EvaluatedResultRounded {
        int32_t distance;
        float weight;
    };

    // the block results share a common prefix layout: the bin indices, followed by the weights, followed by the exact distances for the 
    // non-rounded variants. this lets the kernels fill both variants through the same two pointers, writing the exact distances just past 
    // the end of the weights.
    struct alignas(16) QuadEvaluatedResult {
        std::array<int32_t, 4> distance_bins;
        std::array<float, 4>   weights;
        std::array<float, 4>   distances;
    };

    struct alignas(16) QuadEvaluatedResultRounded {
        std::array<int32_t, 4> distances;
        std::array<float, 4>   weights;
    };

    struct alignas(32) OctoEvaluatedResult {
        std::array<int32_t, 8> distance_bins;
        std::array<float, 8>   weights;
        std::array<float, 8>   distances;
    };

    struct alignas(32) OctoEvaluatedResultRounded {
        std::array<int32_t, 8> distances;
        std::array<float, 8>   weights;
    };

    struct alignas(64) HexaEvaluatedResult {
        std::array<int32_t, 16> distance_bins;
        std::array<float, 16>   weights;
        std::array<float, 16>   distances;
    };

    struct alignas(64) HexaEvaluatedResultRounded {
        std::array<int32_t, 16> distances;
        std::array<float, 16>   weights;
    };

    // assert that it is safe to perform memcpy and reinterpret_cast on these structures
    static_assert(sizeof(EvaluatedResult)            == 12,  "hist::detail::EvaluatedResult is not 12 bytes long");
    static_assert(sizeof(EvaluatedResultRounded)     == 8,   "hist::detail::EvaluatedResultRounded is not 8 bytes long");
    static_assert(sizeof(QuadEvaluatedResult)        == 48,  "hist::detail::QuadEvaluatedResult is not 48 bytes long");
    static_assert(sizeof(QuadEvaluatedResultRounded) == 32,  "hist::detail::QuadEvaluatedResultRounded is not 32 bytes long");
    static_assert(sizeof(OctoEvaluatedResult)        == 96,  "hist::detail::OctoEvaluatedResult is not 96 bytes long");
    static_assert(sizeof(OctoEvaluatedResultRounded) == 64,  "hist::detail::OctoEvaluatedResultRounded is not 64 bytes long");
    static_assert(sizeof(HexaEvaluatedResult)        == 192, "hist::detail::HexaEvaluatedResult is not 192 bytes long");
    static_assert(sizeof(HexaEvaluatedResultRounded) == 128, "hist::detail::HexaEvaluatedResultRounded is not 128 bytes long");

    // ensure our structures are trivially copyable
    static_assert(std::is_trivial_v<EvaluatedResult>,            "hist::detail::EvaluatedResult is not trivial");
    static_assert(std::is_trivial_v<EvaluatedResultRounded>,     "hist::detail::EvaluatedResultRounded is not trivial");
    static_assert(std::is_trivial_v<QuadEvaluatedResult>,        "hist::detail::QuadEvaluatedResult is not trivial");
    static_assert(std::is_trivial_v<QuadEvaluatedResultRounded>, "hist::detail::QuadEvaluatedResultRounded is not trivial");
    static_assert(std::is_trivial_v<OctoEvaluatedResult>,        "hist::detail::OctoEvaluatedResult is not trivial");
    static_assert(std::is_trivial_v<OctoEvaluatedResultRounded>, "hist::detail::OctoEvaluatedResultRounded is not trivial");
    static_assert(std::is_trivial_v<HexaEvaluatedResult>,        "hist::detail::HexaEvaluatedResult is not trivial");
    static_assert(std::is_trivial_v<HexaEvaluatedResultRounded>, "hist::detail::HexaEvaluatedResultRounded is not trivial");

    // check that the structures have a standard memory layout. this is required for the reinterpret_casts.
    static_assert(std::is_standard_layout_v<EvaluatedResult>,            "hist::detail::EvaluatedResult is not standard layout");
    static_assert(std::is_standard_layout_v<EvaluatedResultRounded>,     "hist::detail::EvaluatedResultRounded is not standard layout");
    static_assert(std::is_standard_layout_v<QuadEvaluatedResult>,        "hist::detail::QuadEvaluatedResult is not standard layout");
    static_assert(std::is_standard_layout_v<QuadEvaluatedResultRounded>, "hist::detail::QuadEvaluatedResultRounded is not standard layout");
    static_assert(std::is_standard_layout_v<OctoEvaluatedResult>,        "hist::detail::OctoEvaluatedResult is not standard layout");
    static_assert(std::is_standard_layout_v<OctoEvaluatedResultRounded>, "hist::detail::OctoEvaluatedResultRounded is not standard layout");
    static_assert(std::is_standard_layout_v<HexaEvaluatedResult>,        "hist::detail::HexaEvaluatedResult is not standard layout");
    static_assert(std::is_standard_layout_v<HexaEvaluatedResultRounded>, "hist::detail::HexaEvaluatedResultRounded is not standard layout");

    // the kernels write the exact distances N entries past the start of the weights, so the arrays must be contiguous
    static_assert(offsetof(QuadEvaluatedResult, weights) ==  4*sizeof(float) && offsetof(QuadEvaluatedResult, distances) ==  8*sizeof(float), "hist::detail::QuadEvaluatedResult is not contiguous");
    static_assert(offsetof(OctoEvaluatedResult, weights) ==  8*sizeof(float) && offsetof(OctoEvaluatedResult, distances) == 16*sizeof(float), "hist::detail::OctoEvaluatedResult is not contiguous");
    static_assert(offsetof(HexaEvaluatedResult, weights) == 16*sizeof(float) && offsetof(HexaEvaluatedResult, distances) == 32*sizeof(float), "hist::detail::HexaEvaluatedResult is not contiguous");

    /**
     * @brief A single atom, broadcast against a block of others.
     */
    struct Atom {
        float x = 0, y = 0, z = 0, w = 0;
    };

    /**
     * @brief The first element of a block of atoms, one pointer per component.
     *        The kernels read N consecutive entries from each; the caller guarantees they exist.
     */
    struct Block {
        const float* x = nullptr;
        const float* y = nullptr;
        const float* z = nullptr;
        const float* w = nullptr;
    };

    inline Block advance(Block b, int n) noexcept {return Block{.x=b.x+n, .y=b.y+n, .z=b.z+n, .w=b.w+n};}
}

//#########################################//
//############ IMPLEMENTATION #############//
//#########################################//

// implementation defined in header to support efficient inlining
#if defined AUSAXS_USE_SSE2
    #include <nmmintrin.h>
#endif
#if defined AUSAXS_USE_AVX2 || defined AUSAXS_USE_AVX512
    #include <immintrin.h>
#endif

namespace ausaxs::hist::detail::xyzw {
    //=========================== scalar ===========================//
    /**
     * @brief Evaluate a block of N atoms into the arrays starting at @a bin_out and @a wt_out.
     * 
     * @a W is the width of the result being filled, which is larger than N when a result is assembled from several narrower blocks. 
     * The exact distances are written to the array starting @a W entries past @a wt_out; the default 0 skips them, which is what the 
     * rounded results want.
     */
    template<bool vbw, int N, int W = 0>
    inline void evaluate_N_scalar(Atom self, Block other, int32_t* bin_out, float* wt_out) noexcept {
        const float inv_width = WidthController<vbw>::get_inv_width();
        for (int k = 0; k < N; ++k) {
            float dx = self.x - other.x[k];
            float dy = self.y - other.y[k];
            float dz = self.z - other.z[k];
            float dist = std::sqrt(dx*dx + dy*dy + dz*dz);
            bin_out[k] = static_cast<int32_t>(std::round(inv_width*dist));
            wt_out[k] = self.w*other.w[k];
            if constexpr (W != 0) {(wt_out + W)[k] = dist;}
        }
    }

    //=========================== SSE2 ===========================//
    #if defined AUSAXS_USE_SSE2
        inline void body_4_sse(Atom self, Block other, __m128& dist, __m128& weight) noexcept {
            __m128 dx = _mm_sub_ps(_mm_set_ps1(self.x), _mm_loadu_ps(other.x));
            __m128 dy = _mm_sub_ps(_mm_set_ps1(self.y), _mm_loadu_ps(other.y));
            __m128 dz = _mm_sub_ps(_mm_set_ps1(self.z), _mm_loadu_ps(other.z));
            __m128 d2 = _mm_add_ps(_mm_add_ps(_mm_mul_ps(dx, dx), _mm_mul_ps(dy, dy)), _mm_mul_ps(dz, dz));
            dist = _mm_sqrt_ps(d2);
            weight = _mm_mul_ps(_mm_set_ps1(self.w), _mm_loadu_ps(other.w));
        }

        /// @brief Evaluate a block of 4 atoms. See evaluate_N_scalar for the meaning of the arguments.
        template<bool vbw, int W = 0>
        inline void evaluate_4_sse_into(Atom self, Block other, int32_t* bin_out, float* wt_out) noexcept {
            __m128 dist, weight;
            body_4_sse(self, other, dist, weight);
            _mm_storeu_si128(
                reinterpret_cast<__m128i*>(bin_out),
                _mm_cvtps_epi32(_mm_mul_ps(dist, _mm_set_ps1(WidthController<vbw>::get_inv_width())))
            );
            _mm_storeu_ps(wt_out, weight);
            if constexpr (W != 0) {_mm_storeu_ps(wt_out + W, dist);}
        }
    #endif

    //=========================== AVX2 ===========================//
    #if defined AUSAXS_USE_AVX2
        inline void body_8_avx(Atom self, Block other, __m256& dist, __m256& weight) noexcept {
            __m256 dx = _mm256_sub_ps(_mm256_set1_ps(self.x), _mm256_loadu_ps(other.x));
            __m256 dy = _mm256_sub_ps(_mm256_set1_ps(self.y), _mm256_loadu_ps(other.y));
            __m256 dz = _mm256_sub_ps(_mm256_set1_ps(self.z), _mm256_loadu_ps(other.z));
            __m256 d2 = _mm256_fmadd_ps(dz, dz, _mm256_fmadd_ps(dy, dy, _mm256_mul_ps(dx, dx)));
            dist = _mm256_sqrt_ps(d2);
            weight = _mm256_mul_ps(_mm256_set1_ps(self.w), _mm256_loadu_ps(other.w));
        }

        /// @brief Evaluate a block of 8 atoms. See evaluate_N_scalar for the meaning of the arguments.
        template<bool vbw, int W = 0>
        inline void evaluate_8_avx_into(Atom self, Block other, int32_t* bin_out, float* wt_out) noexcept {
            __m256 dist, weight;
            body_8_avx(self, other, dist, weight);
            _mm256_storeu_si256(
                reinterpret_cast<__m256i*>(bin_out),
                _mm256_cvtps_epi32(_mm256_mul_ps(dist, _mm256_set1_ps(WidthController<vbw>::get_inv_width())))
            );
            _mm256_storeu_ps(wt_out, weight);
            if constexpr (W != 0) {_mm256_storeu_ps(wt_out + W, dist);}
        }
    #endif

    //=========================== AVX512 ===========================//
    #if defined AUSAXS_USE_AVX512
        inline void body_16_avx512(Atom self, Block other, __m512& dist, __m512& weight) noexcept {
            __m512 dx = _mm512_sub_ps(_mm512_set1_ps(self.x), _mm512_loadu_ps(other.x));
            __m512 dy = _mm512_sub_ps(_mm512_set1_ps(self.y), _mm512_loadu_ps(other.y));
            __m512 dz = _mm512_sub_ps(_mm512_set1_ps(self.z), _mm512_loadu_ps(other.z));
            __m512 d2 = _mm512_fmadd_ps(dz, dz, _mm512_fmadd_ps(dy, dy, _mm512_mul_ps(dx, dx)));
            dist = _mm512_sqrt_ps(d2);
            weight = _mm512_mul_ps(_mm512_set1_ps(self.w), _mm512_loadu_ps(other.w));
        }

        /// @brief Evaluate a block of 16 atoms. See evaluate_N_scalar for the meaning of the arguments.
        template<bool vbw, int W = 0>
        inline void evaluate_16_avx512_into(Atom self, Block other, int32_t* bin_out, float* wt_out) noexcept {
            __m512 dist, weight;
            body_16_avx512(self, other, dist, weight);
            _mm512_storeu_si512(
                reinterpret_cast<__m512i*>(bin_out),
                _mm512_cvtps_epi32(_mm512_mul_ps(dist, _mm512_set1_ps(WidthController<vbw>::get_inv_width())))
            );
            _mm512_storeu_ps(wt_out, weight);
            if constexpr (W != 0) {_mm512_storeu_ps(wt_out + W, dist);}
        }
    #endif

    //=========================== dispatch ===========================//
    /**
     * @brief Calculate the distance and combined weight between @a self and a single other atom.
     */
    template<bool vbw>
    inline EvaluatedResult evaluate(Atom self, Block other) noexcept {
        float dx = self.x - other.x[0], dy = self.y - other.y[0], dz = self.z - other.z[0];
        float dist = std::sqrt(dx*dx + dy*dy + dz*dz);
        return EvaluatedResult{
            .distance=dist, 
            .distance_bin=static_cast<int32_t>(std::round(WidthController<vbw>::get_inv_width()*dist)), 
            .weight=self.w*other.w[0]
        };
    }

    /**
     * @brief Calculate the @a binned distance and combined weight between @a self and a single other atom.
     */
    template<bool vbw>
    inline EvaluatedResultRounded evaluate_rounded(Atom self, Block other) noexcept {
        float dx = self.x - other.x[0], dy = self.y - other.y[0], dz = self.z - other.z[0];
        return EvaluatedResultRounded{
            .distance=static_cast<int32_t>(std::round(WidthController<vbw>::get_inv_width()*std::sqrt(dx*dx + dy*dy + dz*dz))),
            .weight=self.w*other.w[0]
        };
    }

    template<bool vbw>
    inline QuadEvaluatedResult evaluate_4(Atom self, Block other) noexcept {
        QuadEvaluatedResult r;
        #if defined AUSAXS_USE_SSE2
            evaluate_4_sse_into<vbw, 4>(self, other, r.distance_bins.data(), r.weights.data());
        #else
            evaluate_N_scalar<vbw, 4, 4>(self, other, r.distance_bins.data(), r.weights.data());
        #endif
        return r;
    }

    template<bool vbw>
    inline QuadEvaluatedResultRounded evaluate_rounded_4(Atom self, Block other) noexcept {
        QuadEvaluatedResultRounded r;
        #if defined AUSAXS_USE_SSE2
            evaluate_4_sse_into<vbw>(self, other, r.distances.data(), r.weights.data());
        #else
            evaluate_N_scalar<vbw, 4>(self, other, r.distances.data(), r.weights.data());
        #endif
        return r;
    }

    template<bool vbw>
    inline OctoEvaluatedResult evaluate_8(Atom self, Block other) noexcept {
        OctoEvaluatedResult r;
        #if defined AUSAXS_USE_AVX2
            evaluate_8_avx_into<vbw, 8>(self, other, r.distance_bins.data(), r.weights.data());
        #elif defined AUSAXS_USE_SSE2
            evaluate_4_sse_into<vbw, 8>(self, other, r.distance_bins.data(), r.weights.data());
            evaluate_4_sse_into<vbw, 8>(self, advance(other, 4), r.distance_bins.data()+4, r.weights.data()+4);
        #else
            evaluate_N_scalar<vbw, 8, 8>(self, other, r.distance_bins.data(), r.weights.data());
        #endif
        return r;
    }

    template<bool vbw>
    inline OctoEvaluatedResultRounded evaluate_rounded_8(Atom self, Block other) noexcept {
        OctoEvaluatedResultRounded r;
        #if defined AUSAXS_USE_AVX2
            evaluate_8_avx_into<vbw>(self, other, r.distances.data(), r.weights.data());
        #elif defined AUSAXS_USE_SSE2
            evaluate_4_sse_into<vbw>(self, other, r.distances.data(), r.weights.data());
            evaluate_4_sse_into<vbw>(self, advance(other, 4), r.distances.data()+4, r.weights.data()+4);
        #else
            evaluate_N_scalar<vbw, 8>(self, other, r.distances.data(), r.weights.data());
        #endif
        return r;
    }

    template<bool vbw>
    inline HexaEvaluatedResult evaluate_16(Atom self, Block other) noexcept {
        HexaEvaluatedResult r;
        #if defined AUSAXS_USE_AVX512
            evaluate_16_avx512_into<vbw, 16>(self, other, r.distance_bins.data(), r.weights.data());
        #elif defined AUSAXS_USE_AVX2
            evaluate_8_avx_into<vbw, 16>(self, other, r.distance_bins.data(), r.weights.data());
            evaluate_8_avx_into<vbw, 16>(self, advance(other, 8), r.distance_bins.data()+8, r.weights.data()+8);
        #elif defined AUSAXS_USE_SSE2
            for (int b = 0; b < 4; ++b) {
                evaluate_4_sse_into<vbw, 16>(self, advance(other, 4*b), r.distance_bins.data()+4*b, r.weights.data()+4*b);
            }
        #else
            evaluate_N_scalar<vbw, 16, 16>(self, other, r.distance_bins.data(), r.weights.data());
        #endif
        return r;
    }

    template<bool vbw>
    inline HexaEvaluatedResultRounded evaluate_rounded_16(Atom self, Block other) noexcept {
        HexaEvaluatedResultRounded r;
        #if defined AUSAXS_USE_AVX512
            evaluate_16_avx512_into<vbw>(self, other, r.distances.data(), r.weights.data());
        #elif defined AUSAXS_USE_AVX2
            evaluate_8_avx_into<vbw>(self, other, r.distances.data(), r.weights.data());
            evaluate_8_avx_into<vbw>(self, advance(other, 8), r.distances.data()+8, r.weights.data()+8);
        #elif defined AUSAXS_USE_SSE2
            for (int b = 0; b < 4; ++b) {
                evaluate_4_sse_into<vbw>(self, advance(other, 4*b), r.distances.data()+4*b, r.weights.data()+4*b);
            }
        #else
            evaluate_N_scalar<vbw, 16>(self, other, r.distances.data(), r.weights.data());
        #endif
        return r;
    }
}
