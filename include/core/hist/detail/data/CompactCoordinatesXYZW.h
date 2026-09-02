// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

/**
 * @brief This file contains custom vector instructions for efficient scattering calculations.
 *
 * The implementation is specialized for generic systems defined by separate x, y, z and w arrays.
 * This is useful for the Simple excluded volume model, where all atoms have the same form factor
 * type but different weights. Similarly, it is useful for SANS calculations, where the lack of
 * q-dependence means the form factors may be encoded as simple weights.
 *
 * For more complex X-ray calculations with different atomic species, the CompactCoordinatesXYZFF
 * implementation may be more useful, since it is specialized for systems where the fourth
 * component is an int32 form factor index.
 *
 * The coordinates are passed as one pointer per component rather than as interleaved [x, y, z, w]
 * tuples. Interleaved storage forces every kernel to transpose a block before it can compute the
 * squared distance, which cost 8-10 shuffle-port operations per block; with separate arrays there
 * is nothing to transpose and the kernels reduce to loads, subtractions and multiply-adds. That
 * was worth 24-45% of the kernel depending on dispatch level, and it also lets the compiler
 * auto-vectorise the fallback path competitively, which matters on targets with no hand-written
 * SIMD at all.
 */

#pragma once

#include <hist/detail/data/IntrinsicMacros.h>
#include <hist/detail/data/WidthControllers.h>
#include <hist/detail/data/IntrinsicHelpers.h>
#include <constants/Constants.h>
#include <settings/Flags.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <type_traits>

namespace ausaxs::hist::detail::xyzw {
    struct EvaluatedResult {
        float distance;         // The raw distance
        int32_t distance_bin;   // The distance bin index
        float weight;           // The combined weight
    };

    struct EvaluatedResultRounded {
        int32_t distance;       // The distance bin
        float weight;           // The combined weight
    };

    struct alignas(16) QuadEvaluatedResult {
        std::array<float, 4> distances;       // The raw distances (for weighted bin center calculation)
        std::array<int32_t, 4> distance_bins; // The distance bin indices (for array indexing)
        std::array<float, 4> weights;         // The combined weight
    };

    struct alignas(16) QuadEvaluatedResultRounded {
        std::array<int32_t, 4> distances;   // The distance bin
        std::array<float, 4> weights;       // The combined weight
    };

    struct alignas(32) OctoEvaluatedResult {
        std::array<float, 8> distances;       // The raw distances (for weighted bin center calculation)
        std::array<int32_t, 8> distance_bins; // The distance bin indices (for array indexing)
        std::array<float, 8> weights;         // The combined weight
    };

    struct alignas(32) OctoEvaluatedResultRounded {
        std::array<int32_t, 8> distances;   // The distance bin
        std::array<float, 8> weights;       // The combined weight
    };

    // 64-byte aligned, not 32: each field is a whole 512-bit register's worth, and at 32-byte
    // alignment the second field lands on a 32-mod-64 address. The kernels below use the
    // unaligned store intrinsics so neither this nor the sub-block offsets can fault, but the
    // alignment is what makes those stores single-cache-line in practice. Do not lower it.
    struct alignas(64) HexaEvaluatedResult {
        std::array<float, 16> distances;
        std::array<int32_t, 16> distance_bins;
        std::array<float, 16> weights;
    };

    struct alignas(64) HexaEvaluatedResultRounded {
        std::array<int32_t, 16> distances;
        std::array<float, 16> weights;
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

    inline Block advance(Block b, int n) noexcept {return Block{b.x + n, b.y + n, b.z + n, b.w + n};}
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
    // Deliberately free of intrinsics: with the coordinates in separate arrays the compiler
    // vectorises this on its own, so this is the portable path rather than a slow fallback.
    template<bool vbw, int N, typename Result>
    inline void evaluate_N_scalar(Atom self, Block other, Result& out) noexcept {
        const float inv_width = WidthController<vbw>::get_inv_width();
        for (int k = 0; k < N; ++k) {
            float dx = self.x - other.x[k];
            float dy = self.y - other.y[k];
            float dz = self.z - other.z[k];
            float dist = std::sqrt(dx*dx + dy*dy + dz*dz);
            out.distances[k] = dist;
            out.distance_bins[k] = static_cast<int32_t>(std::round(inv_width*dist));
            out.weights[k] = self.w*other.w[k];
        }
    }

    template<bool vbw, int N, typename Result>
    inline void evaluate_rounded_N_scalar(Atom self, Block other, Result& out) noexcept {
        const float inv_width = WidthController<vbw>::get_inv_width();
        for (int k = 0; k < N; ++k) {
            float dx = self.x - other.x[k];
            float dy = self.y - other.y[k];
            float dz = self.z - other.z[k];
            out.distances[k] = static_cast<int32_t>(std::round(inv_width*std::sqrt(dx*dx + dy*dy + dz*dz)));
            out.weights[k] = self.w*other.w[k];
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

        template<bool vbw>
        inline void evaluate_4_sse_into(Atom self, Block other, float* dist_out, int32_t* bin_out, float* wt_out) noexcept {
            __m128 dist, weight;
            body_4_sse(self, other, dist, weight);
            _mm_storeu_ps(dist_out, dist);
            _mm_storeu_si128(reinterpret_cast<__m128i*>(bin_out),
                _mm_cvtps_epi32(_mm_mul_ps(dist, _mm_set_ps1(WidthController<vbw>::get_inv_width()))));
            _mm_storeu_ps(wt_out, weight);
        }

        template<bool vbw>
        inline void evaluate_rounded_4_sse_into(Atom self, Block other, int32_t* dist_out, float* wt_out) noexcept {
            __m128 dist, weight;
            body_4_sse(self, other, dist, weight);
            _mm_storeu_si128(reinterpret_cast<__m128i*>(dist_out),
                _mm_cvtps_epi32(_mm_mul_ps(dist, _mm_set_ps1(WidthController<vbw>::get_inv_width()))));
            _mm_storeu_ps(wt_out, weight);
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

        template<bool vbw>
        inline void evaluate_8_avx_into(Atom self, Block other, float* dist_out, int32_t* bin_out, float* wt_out) noexcept {
            __m256 dist, weight;
            body_8_avx(self, other, dist, weight);
            _mm256_storeu_ps(dist_out, dist);
            _mm256_storeu_si256(reinterpret_cast<__m256i*>(bin_out),
                _mm256_cvtps_epi32(_mm256_mul_ps(dist, _mm256_set1_ps(WidthController<vbw>::get_inv_width()))));
            _mm256_storeu_ps(wt_out, weight);
        }

        template<bool vbw>
        inline void evaluate_rounded_8_avx_into(Atom self, Block other, int32_t* dist_out, float* wt_out) noexcept {
            __m256 dist, weight;
            body_8_avx(self, other, dist, weight);
            _mm256_storeu_si256(reinterpret_cast<__m256i*>(dist_out),
                _mm256_cvtps_epi32(_mm256_mul_ps(dist, _mm256_set1_ps(WidthController<vbw>::get_inv_width()))));
            _mm256_storeu_ps(wt_out, weight);
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

        template<bool vbw>
        inline void evaluate_16_avx512_into(Atom self, Block other, float* dist_out, int32_t* bin_out, float* wt_out) noexcept {
            __m512 dist, weight;
            body_16_avx512(self, other, dist, weight);
            _mm512_storeu_ps(dist_out, dist);
            _mm512_storeu_si512(reinterpret_cast<__m512i*>(bin_out),
                _mm512_cvtps_epi32(_mm512_mul_ps(dist, _mm512_set1_ps(WidthController<vbw>::get_inv_width()))));
            _mm512_storeu_ps(wt_out, weight);
        }

        template<bool vbw>
        inline void evaluate_rounded_16_avx512_into(Atom self, Block other, int32_t* dist_out, float* wt_out) noexcept {
            __m512 dist, weight;
            body_16_avx512(self, other, dist, weight);
            _mm512_storeu_si512(reinterpret_cast<__m512i*>(dist_out),
                _mm512_cvtps_epi32(_mm512_mul_ps(dist, _mm512_set1_ps(WidthController<vbw>::get_inv_width()))));
            _mm512_storeu_ps(wt_out, weight);
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
        return EvaluatedResult{dist, static_cast<int32_t>(std::round(WidthController<vbw>::get_inv_width()*dist)), self.w*other.w[0]};
    }

    /**
     * @brief Calculate the @a binned distance and combined weight between @a self and a single other atom.
     */
    template<bool vbw>
    inline EvaluatedResultRounded evaluate_rounded(Atom self, Block other) noexcept {
        float dx = self.x - other.x[0], dy = self.y - other.y[0], dz = self.z - other.z[0];
        return EvaluatedResultRounded{
            static_cast<int32_t>(std::round(WidthController<vbw>::get_inv_width()*std::sqrt(dx*dx + dy*dy + dz*dz))),
            self.w*other.w[0]
        };
    }

    template<bool vbw>
    inline QuadEvaluatedResult evaluate_4(Atom self, Block other) noexcept {
        QuadEvaluatedResult r;
        #if defined AUSAXS_USE_SSE2
            evaluate_4_sse_into<vbw>(self, other, r.distances.data(), r.distance_bins.data(), r.weights.data());
        #else
            evaluate_N_scalar<vbw, 4>(self, other, r);
        #endif
        return r;
    }

    template<bool vbw>
    inline QuadEvaluatedResultRounded evaluate_rounded_4(Atom self, Block other) noexcept {
        QuadEvaluatedResultRounded r;
        #if defined AUSAXS_USE_SSE2
            evaluate_rounded_4_sse_into<vbw>(self, other, r.distances.data(), r.weights.data());
        #else
            evaluate_rounded_N_scalar<vbw, 4>(self, other, r);
        #endif
        return r;
    }

    template<bool vbw>
    inline OctoEvaluatedResult evaluate_8(Atom self, Block other) noexcept {
        OctoEvaluatedResult r;
        #if defined AUSAXS_USE_AVX2
            evaluate_8_avx_into<vbw>(self, other, r.distances.data(), r.distance_bins.data(), r.weights.data());
        #elif defined AUSAXS_USE_SSE2
            evaluate_4_sse_into<vbw>(self, other, r.distances.data(), r.distance_bins.data(), r.weights.data());
            evaluate_4_sse_into<vbw>(self, advance(other, 4), r.distances.data()+4, r.distance_bins.data()+4, r.weights.data()+4);
        #else
            evaluate_N_scalar<vbw, 8>(self, other, r);
        #endif
        return r;
    }

    template<bool vbw>
    inline OctoEvaluatedResultRounded evaluate_rounded_8(Atom self, Block other) noexcept {
        OctoEvaluatedResultRounded r;
        #if defined AUSAXS_USE_AVX2
            evaluate_rounded_8_avx_into<vbw>(self, other, r.distances.data(), r.weights.data());
        #elif defined AUSAXS_USE_SSE2
            evaluate_rounded_4_sse_into<vbw>(self, other, r.distances.data(), r.weights.data());
            evaluate_rounded_4_sse_into<vbw>(self, advance(other, 4), r.distances.data()+4, r.weights.data()+4);
        #else
            evaluate_rounded_N_scalar<vbw, 8>(self, other, r);
        #endif
        return r;
    }

    template<bool vbw>
    inline HexaEvaluatedResult evaluate_16(Atom self, Block other) noexcept {
        HexaEvaluatedResult r;
        #if defined AUSAXS_USE_AVX512
            evaluate_16_avx512_into<vbw>(self, other, r.distances.data(), r.distance_bins.data(), r.weights.data());
        #elif defined AUSAXS_USE_AVX2
            evaluate_8_avx_into<vbw>(self, other, r.distances.data(), r.distance_bins.data(), r.weights.data());
            evaluate_8_avx_into<vbw>(self, advance(other, 8), r.distances.data()+8, r.distance_bins.data()+8, r.weights.data()+8);
        #elif defined AUSAXS_USE_SSE2
            for (int b = 0; b < 4; ++b) {
                evaluate_4_sse_into<vbw>(self, advance(other, 4*b), r.distances.data()+4*b, r.distance_bins.data()+4*b, r.weights.data()+4*b);
            }
        #else
            evaluate_N_scalar<vbw, 16>(self, other, r);
        #endif
        return r;
    }

    template<bool vbw>
    inline HexaEvaluatedResultRounded evaluate_rounded_16(Atom self, Block other) noexcept {
        HexaEvaluatedResultRounded r;
        #if defined AUSAXS_USE_AVX512
            evaluate_rounded_16_avx512_into<vbw>(self, other, r.distances.data(), r.weights.data());
        #elif defined AUSAXS_USE_AVX2
            evaluate_rounded_8_avx_into<vbw>(self, other, r.distances.data(), r.weights.data());
            evaluate_rounded_8_avx_into<vbw>(self, advance(other, 8), r.distances.data()+8, r.weights.data()+8);
        #elif defined AUSAXS_USE_SSE2
            for (int b = 0; b < 4; ++b) {
                evaluate_rounded_4_sse_into<vbw>(self, advance(other, 4*b), r.distances.data()+4*b, r.weights.data()+4*b);
            }
        #else
            evaluate_rounded_N_scalar<vbw, 16>(self, other, r);
        #endif
        return r;
    }
}
