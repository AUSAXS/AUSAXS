// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

/**
 * @brief This file contains custom vector instructions for efficient scattering calculations.
 * 
 * The implementation is specialized for generic systems defined by a collection of [x, y, z, w] float vectors.
 * This is useful for the Simple excluded volume model, where all atoms have the same form factor type but different weights.
 * Similarly, it is useful for SANS calculations, where the lack of q-dependence means the form factors may be encoded as simple weights.
 * 
 * For more complex X-ray calculations with different atomic species, the CompactCoordinatesXYZWFF implementation may be more useful,
 * since it is specialized for systems defined by a collection of [x: float, y: float, z: float, ff_type: int32] vectors.
 */

 #pragma once

#include <hist/detail/data/IntrinsicMacros.h>
#include <hist/detail/data/WidthControllers.h>
#include <hist/detail/data/IntrinsicHelpers.h>
#include <math/Vector3.h>
#include <constants/Constants.h>
#include <settings/Flags.h>

#include <array>
#include <cstdint>
#include <span>

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

    struct QuadEvaluatedResult {
        std::array<float, 4> distances;       // The raw distances (for weighted bin center calculation)
        std::array<int32_t, 4> distance_bins; // The distance bin indices (for array indexing)
        std::array<float, 4> weights;         // The combined weight
    };

    struct QuadEvaluatedResultRounded {
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
    static_assert(std::is_standard_layout_v<EvaluatedResult>,            "hist::detail::EvaluatedResult is not trivial");
    static_assert(std::is_standard_layout_v<EvaluatedResultRounded>,     "hist::detail::EvaluatedResultRounded is not trivial");
    static_assert(std::is_standard_layout_v<QuadEvaluatedResult>,        "hist::detail::QuadEvaluatedResult is not trivial");
    static_assert(std::is_standard_layout_v<QuadEvaluatedResultRounded>, "hist::detail::QuadEvaluatedResultRounded is not trivial");
    static_assert(std::is_standard_layout_v<OctoEvaluatedResult>,        "hist::detail::OctoEvaluatedResult is not trivial");
    static_assert(std::is_standard_layout_v<OctoEvaluatedResultRounded>, "hist::detail::OctoEvaluatedResultRounded is not trivial");
    static_assert(std::is_standard_layout_v<HexaEvaluatedResult>,        "hist::detail::HexaEvaluatedResult is not trivial");
    static_assert(std::is_standard_layout_v<HexaEvaluatedResultRounded>, "hist::detail::HexaEvaluatedResultRounded is not trivial");
}

namespace ausaxs::hist::detail {
    template<bool variable_bin_width>
    class CompactCoordinatesXYZW : public WidthController<variable_bin_width> {
        public:
            using WidthController<variable_bin_width>::get_inv_width;
            CompactCoordinatesXYZW() noexcept = default;
            CompactCoordinatesXYZW(const CompactCoordinatesXYZW& other) noexcept = default;
            CompactCoordinatesXYZW(CompactCoordinatesXYZW&& other) noexcept = default;
            CompactCoordinatesXYZW& operator=(const CompactCoordinatesXYZW& other) noexcept = default;
            CompactCoordinatesXYZW& operator=(CompactCoordinatesXYZW&& other) noexcept = default;

            template<numeric T, numeric V>
            CompactCoordinatesXYZW(const Vector3<T>& v, V w) noexcept : value{.pos={static_cast<float>(v.x()), static_cast<float>(v.y()), static_cast<float>(v.z())}, .w=static_cast<float>(w)} {}
            CompactCoordinatesXYZW(const Vector3<float>& v, float w) noexcept : value{.pos=v, .w=w} {}

            xyzw::EvaluatedResultRounded evaluate_rounded(const CompactCoordinatesXYZW& other) const noexcept;
            xyzw::EvaluatedResult evaluate(const CompactCoordinatesXYZW& other) const noexcept;

            xyzw::QuadEvaluatedResultRounded evaluate_rounded_4(std::span<const CompactCoordinatesXYZW, 4> others) const noexcept;
            xyzw::QuadEvaluatedResult evaluate_4(std::span<const CompactCoordinatesXYZW, 4> others) const noexcept;

            xyzw::OctoEvaluatedResultRounded evaluate_rounded_8(std::span<const CompactCoordinatesXYZW, 8> others) const noexcept;
            xyzw::OctoEvaluatedResult evaluate_8(std::span<const CompactCoordinatesXYZW, 8> others) const noexcept;

            xyzw::HexaEvaluatedResultRounded evaluate_rounded_16(std::span<const CompactCoordinatesXYZW, 16> others) const noexcept;
            xyzw::HexaEvaluatedResult evaluate_16(std::span<const CompactCoordinatesXYZW, 16> others) const noexcept;

            union {
                struct {Vector3<float> pos; float w;} value;
                std::array<float, 4> data;
            };

        protected:
            xyzw::EvaluatedResultRounded evaluate_rounded_scalar(const CompactCoordinatesXYZW& other) const noexcept;
            xyzw::EvaluatedResult evaluate_scalar(const CompactCoordinatesXYZW& other) const noexcept;

            xyzw::QuadEvaluatedResultRounded evaluate_rounded_4_scalar(std::span<const CompactCoordinatesXYZW, 4> others) const noexcept;
            xyzw::QuadEvaluatedResult evaluate_4_scalar(std::span<const CompactCoordinatesXYZW, 4> others) const noexcept;

            #if defined AUSAXS_USE_SSE2
                xyzw::QuadEvaluatedResultRounded evaluate_rounded_4_sse(std::span<const CompactCoordinatesXYZW, 4> others) const noexcept;
                xyzw::QuadEvaluatedResult evaluate_4_sse(std::span<const CompactCoordinatesXYZW, 4> others) const noexcept;
                void evaluate_rounded_4_sse_into(std::span<const CompactCoordinatesXYZW, 4> others, int32_t* dist_out, float* wt_out) const noexcept;
                void evaluate_4_sse_into(std::span<const CompactCoordinatesXYZW, 4> others, float* dist_out, int32_t* bin_out, float* wt_out) const noexcept;
            #endif

            #if defined AUSAXS_USE_AVX2
                xyzw::OctoEvaluatedResultRounded evaluate_rounded_8_avx(std::span<const CompactCoordinatesXYZW, 8> others) const noexcept;
                xyzw::OctoEvaluatedResult evaluate_8_avx(std::span<const CompactCoordinatesXYZW, 8> others) const noexcept;
                void evaluate_rounded_8_avx_into(std::span<const CompactCoordinatesXYZW, 8> others, int32_t* dist_out, float* wt_out) const noexcept;
                void evaluate_8_avx_into(std::span<const CompactCoordinatesXYZW, 8> others, float* dist_out, int32_t* bin_out, float* wt_out) const noexcept;
            #endif

            #if defined AUSAXS_USE_AVX512
                xyzw::HexaEvaluatedResultRounded evaluate_rounded_16_avx512(std::span<const CompactCoordinatesXYZW, 16> others) const noexcept;
                xyzw::HexaEvaluatedResult evaluate_16_avx512(std::span<const CompactCoordinatesXYZW, 16> others) const noexcept;
            #endif
    };
    static_assert(sizeof(CompactCoordinatesXYZW<true>) == 16,              "CompactCoordinatesXYZW is not 16 bytes. This is required for aligning SIMD instructions.");
    static_assert(std::is_trivial_v<CompactCoordinatesXYZW<true>>,         "CompactCoordinatesXYZW is not trivial");
    static_assert(std::is_standard_layout_v<CompactCoordinatesXYZW<true>>, "CompactCoordinatesXYZW is not standard layout");
    static_assert(supports_nothrow_move_v<CompactCoordinatesXYZW<true>>,   "CompactCoordinatesXYZW should support nothrow move semantics.");
    static_assert(sizeof(CompactCoordinatesXYZW<false>) == 16,             "CompactCoordinatesXYZW is not 16 bytes. This is required for aligning SIMD instructions.");
    static_assert(std::is_trivial_v<CompactCoordinatesXYZW<false>>,        "CompactCoordinatesXYZW is not trivial");
    static_assert(std::is_standard_layout_v<CompactCoordinatesXYZW<false>>,"CompactCoordinatesXYZW is not standard layout");
    static_assert(supports_nothrow_move_v<CompactCoordinatesXYZW<false>>,  "CompactCoordinatesXYZW should support nothrow move semantics.");
}

//#########################################//
//############ IMPLEMENTATION #############//
//#########################################//

// implementation defined in header to support efficient inlining

template<bool vbw>
inline ausaxs::hist::detail::xyzw::EvaluatedResult ausaxs::hist::detail::CompactCoordinatesXYZW<vbw>::evaluate(const CompactCoordinatesXYZW& other) const noexcept {
    return evaluate_scalar(other);
}

template<bool vbw>
inline ausaxs::hist::detail::xyzw::EvaluatedResultRounded ausaxs::hist::detail::CompactCoordinatesXYZW<vbw>::evaluate_rounded(const CompactCoordinatesXYZW& other) const noexcept {
    return evaluate_rounded_scalar(other);
}

template<bool vbw>
inline ausaxs::hist::detail::xyzw::QuadEvaluatedResult ausaxs::hist::detail::CompactCoordinatesXYZW<vbw>::evaluate_4(
    std::span<const CompactCoordinatesXYZW, 4> others
) const noexcept {
    #if defined AUSAXS_USE_SSE2
        return evaluate_4_sse(others);
    #else
        return evaluate_4_scalar(others);
    #endif
}

template<bool vbw>
inline ausaxs::hist::detail::xyzw::QuadEvaluatedResultRounded ausaxs::hist::detail::CompactCoordinatesXYZW<vbw>::evaluate_rounded_4(
    std::span<const CompactCoordinatesXYZW, 4> others
) const noexcept {
    #if defined AUSAXS_USE_SSE2
        return evaluate_rounded_4_sse(others);
    #else
        return evaluate_rounded_4_scalar(others);
    #endif
}

template<bool vbw>
inline ausaxs::hist::detail::xyzw::OctoEvaluatedResult ausaxs::hist::detail::CompactCoordinatesXYZW<vbw>::evaluate_8(
    std::span<const CompactCoordinatesXYZW, 8> others
) const noexcept {
    #if defined AUSAXS_USE_AVX2
        return evaluate_8_avx(others);
    #elif defined AUSAXS_USE_SSE2
        xyzw::OctoEvaluatedResult result;
        evaluate_4_sse_into(others.template first<4>(), result.distances.data(), result.distance_bins.data(), result.weights.data());
        evaluate_4_sse_into(others.template last<4>(), result.distances.data()+4, result.distance_bins.data()+4, result.weights.data()+4);
        return result;
    #else
        return xyzw::OctoEvaluatedResult(
            evaluate_4_scalar(others.template first<4>()),
            evaluate_4_scalar(others.template last<4>())
        );
    #endif
}

template<bool vbw>
inline ausaxs::hist::detail::xyzw::OctoEvaluatedResultRounded ausaxs::hist::detail::CompactCoordinatesXYZW<vbw>::evaluate_rounded_8(
    std::span<const CompactCoordinatesXYZW, 8> others
) const noexcept {
    #if defined AUSAXS_USE_AVX2
        return evaluate_rounded_8_avx(others);
    #elif defined AUSAXS_USE_SSE2
        xyzw::OctoEvaluatedResultRounded result;
        evaluate_rounded_4_sse_into(others.template first<4>(), result.distances.data(), result.weights.data());
        evaluate_rounded_4_sse_into(others.template last<4>(), result.distances.data()+4, result.weights.data()+4);
        return result;
    #else
        return xyzw::OctoEvaluatedResultRounded(
            evaluate_rounded_4_scalar(others.template first<4>()),
            evaluate_rounded_4_scalar(others.template last<4>())
        );
    #endif
}

template<bool vbw>
inline ausaxs::hist::detail::xyzw::HexaEvaluatedResult ausaxs::hist::detail::CompactCoordinatesXYZW<vbw>::evaluate_16(
    std::span<const CompactCoordinatesXYZW, 16> others
) const noexcept {
    #if defined AUSAXS_USE_AVX512
        return evaluate_16_avx512(others);
    #elif defined AUSAXS_USE_AVX2
        xyzw::HexaEvaluatedResult result;
        evaluate_8_avx_into(others.template first<8>(), result.distances.data(), result.distance_bins.data(), result.weights.data());
        evaluate_8_avx_into(others.template last<8>(), result.distances.data()+8, result.distance_bins.data()+8, result.weights.data()+8);
        return result;
    #elif defined AUSAXS_USE_SSE2
        xyzw::HexaEvaluatedResult result;
        evaluate_4_sse_into(others.template first<4>(), result.distances.data(), result.distance_bins.data(), result.weights.data());
        evaluate_4_sse_into(others.template subspan<4,4>(), result.distances.data()+4, result.distance_bins.data()+4, result.weights.data()+4);
        evaluate_4_sse_into(others.template subspan<8,4>(), result.distances.data()+8, result.distance_bins.data()+8, result.weights.data()+8);
        evaluate_4_sse_into(others.template last<4>(), result.distances.data()+12, result.distance_bins.data()+12, result.weights.data()+12);
        return result;
    #else
        return xyzw::HexaEvaluatedResult(evaluate_8(others.template first<8>()), evaluate_8(others.template last<8>()));
    #endif
}

template<bool vbw>
inline ausaxs::hist::detail::xyzw::HexaEvaluatedResultRounded ausaxs::hist::detail::CompactCoordinatesXYZW<vbw>::evaluate_rounded_16(
    std::span<const CompactCoordinatesXYZW, 16> others
) const noexcept {
    #if defined AUSAXS_USE_AVX512
        return evaluate_rounded_16_avx512(others);
    #elif defined AUSAXS_USE_AVX2
        xyzw::HexaEvaluatedResultRounded result;
        evaluate_rounded_8_avx_into(others.template first<8>(), result.distances.data(), result.weights.data());
        evaluate_rounded_8_avx_into(others.template last<8>(), result.distances.data()+8, result.weights.data()+8);
        return result;
    #elif defined AUSAXS_USE_SSE2
        xyzw::HexaEvaluatedResultRounded result;
        evaluate_rounded_4_sse_into(others.template first<4>(), result.distances.data(), result.weights.data());
        evaluate_rounded_4_sse_into(others.template subspan<4,4>(), result.distances.data()+4, result.weights.data()+4);
        evaluate_rounded_4_sse_into(others.template subspan<8,4>(), result.distances.data()+8, result.weights.data()+8);
        evaluate_rounded_4_sse_into(others.template last<4>(), result.distances.data()+12, result.weights.data()+12);
        return result;
    #else
        return xyzw::HexaEvaluatedResultRounded(evaluate_rounded_8(others.template first<8>()), evaluate_rounded_8(others.template last<8>()));
    #endif
}

template<bool vbw>
inline ausaxs::hist::detail::xyzw::EvaluatedResult ausaxs::hist::detail::CompactCoordinatesXYZW<vbw>::evaluate_scalar(const CompactCoordinatesXYZW& other) const noexcept {
    float dist = std::sqrt(squared_dot_product(this->data.data(), other.data.data()));
    int32_t dist_bin = std::round(get_inv_width() * dist);
    return xyzw::EvaluatedResult(dist, dist_bin, value.w*other.value.w);
}

template<bool vbw>
inline ausaxs::hist::detail::xyzw::EvaluatedResultRounded ausaxs::hist::detail::CompactCoordinatesXYZW<vbw>::evaluate_rounded_scalar(
    const CompactCoordinatesXYZW& other
) const noexcept {
    int32_t dist = std::round(get_inv_width()*std::sqrt(squared_dot_product(this->data.data(), other.data.data())));
    return xyzw::EvaluatedResultRounded(dist, value.w*other.value.w);
}

template<bool vbw>
inline ausaxs::hist::detail::xyzw::QuadEvaluatedResult ausaxs::hist::detail::CompactCoordinatesXYZW<vbw>::evaluate_4_scalar(
    std::span<const CompactCoordinatesXYZW, 4> others
) const noexcept {
    const auto& v1 = others[0]; const auto& v2 = others[1]; const auto& v3 = others[2]; const auto& v4 = others[3];
    float dx1 = std::sqrt(squared_dot_product(this->data.data(), v1.data.data()));
    float dx2 = std::sqrt(squared_dot_product(this->data.data(), v2.data.data()));
    float dx3 = std::sqrt(squared_dot_product(this->data.data(), v3.data.data()));
    float dx4 = std::sqrt(squared_dot_product(this->data.data(), v4.data.data()));
    float inv_width = get_inv_width();
    return xyzw::QuadEvaluatedResult(
        std::array<float, 4>{dx1, dx2, dx3, dx4},
        std::array<int32_t, 4>{
            static_cast<int32_t>(std::round(inv_width * dx1)),
            static_cast<int32_t>(std::round(inv_width * dx2)),
            static_cast<int32_t>(std::round(inv_width * dx3)),
            static_cast<int32_t>(std::round(inv_width * dx4))
        },
        std::array<float, 4>{value.w*v1.value.w, value.w*v2.value.w, value.w*v3.value.w, value.w*v4.value.w}
    );
}

template<bool vbw>
inline ausaxs::hist::detail::xyzw::QuadEvaluatedResultRounded ausaxs::hist::detail::CompactCoordinatesXYZW<vbw>::evaluate_rounded_4_scalar(
    std::span<const CompactCoordinatesXYZW, 4> others
) const noexcept {
    const auto& v1 = others[0]; const auto& v2 = others[1]; const auto& v3 = others[2]; const auto& v4 = others[3];
    int32_t dx1 = std::round(get_inv_width()*std::sqrt(squared_dot_product(this->data.data(), v1.data.data())));
    int32_t dx2 = std::round(get_inv_width()*std::sqrt(squared_dot_product(this->data.data(), v2.data.data())));
    int32_t dx3 = std::round(get_inv_width()*std::sqrt(squared_dot_product(this->data.data(), v3.data.data())));
    int32_t dx4 = std::round(get_inv_width()*std::sqrt(squared_dot_product(this->data.data(), v4.data.data())));
    return xyzw::QuadEvaluatedResultRounded(
        std::array<int32_t, 4>{dx1, dx2, dx3, dx4},
        std::array<float, 4>{value.w*v1.value.w, value.w*v2.value.w, value.w*v3.value.w, value.w*v4.value.w}
    );
}

#if defined AUSAXS_USE_SSE2
    #include <nmmintrin.h>
    template<bool vbw>
    inline void ausaxs::hist::detail::CompactCoordinatesXYZW<vbw>::evaluate_4_sse_into(
        std::span<const CompactCoordinatesXYZW, 4> others,
        float* dist_out, int32_t* bin_out, float* wt_out
    ) const noexcept {
        const auto& v1 = others[0]; const auto& v2 = others[1]; const auto& v3 = others[2]; const auto& v4 = others[3];
        __m128 sv = _mm_load_ps(this->data.data());
        __m128 d1 = _mm_sub_ps(sv, _mm_load_ps(v1.data.data()));
        __m128 d2 = _mm_sub_ps(sv, _mm_load_ps(v2.data.data()));
        __m128 d3 = _mm_sub_ps(sv, _mm_load_ps(v3.data.data()));
        __m128 d4 = _mm_sub_ps(sv, _mm_load_ps(v4.data.data()));
        d1 = _mm_mul_ps(d1, d1);
        d2 = _mm_mul_ps(d2, d2);
        d3 = _mm_mul_ps(d3, d3);
        d4 = _mm_mul_ps(d4, d4);

        // 4x4 transpose: converts per-atom [dx²,dy²,dz²,dw²] to per-component rows, discarding w
        __m128 t0 = _mm_unpacklo_ps(d1, d2);   // [dx1²,dx2²,dy1²,dy2²]
        __m128 t1 = _mm_unpackhi_ps(d1, d2);   // [dz1²,dz2²,dw1²,dw2²]
        __m128 t2 = _mm_unpacklo_ps(d3, d4);   // [dx3²,dx4²,dy3²,dy4²]
        __m128 t3 = _mm_unpackhi_ps(d3, d4);   // [dz3²,dz4²,dw3²,dw4²]
        __m128 row_x = _mm_movelh_ps(t0, t2);  // [dx1²,dx2²,dx3²,dx4²]
        __m128 row_y = _mm_movehl_ps(t2, t0);  // [dy1²,dy2²,dy3²,dy4²]
        __m128 row_z = _mm_movelh_ps(t1, t3);  // [dz1²,dz2²,dz3²,dz4²]

        __m128 dist2 = _mm_add_ps(_mm_add_ps(row_x, row_y), row_z);
        __m128 dist_sqrt = _mm_sqrt_ps(dist2);
        __m128 dist_binf = _mm_mul_ps(dist_sqrt, _mm_set_ps1(get_inv_width()));
        __m128i dist_bin = _mm_cvtps_epi32(dist_binf);

        __m128 w1 = _mm_set_ps1(value.w);
        __m128 w2 = _mm_set_ps(v4.value.w, v3.value.w, v2.value.w, v1.value.w);
        __m128 weights = _mm_mul_ps(w1, w2);

        _mm_store_ps(dist_out, dist_sqrt);
        _mm_store_si128(reinterpret_cast<__m128i*>(bin_out), dist_bin);
        _mm_store_ps(wt_out, weights);
    }

    template<bool vbw>
    inline ausaxs::hist::detail::xyzw::QuadEvaluatedResult ausaxs::hist::detail::CompactCoordinatesXYZW<vbw>::evaluate_4_sse(
        std::span<const CompactCoordinatesXYZW, 4> others
    ) const noexcept {
        xyzw::QuadEvaluatedResult result;
        evaluate_4_sse_into(others, result.distances.data(), result.distance_bins.data(), result.weights.data());
        return result;
    }

    template<bool vbw>
    inline void ausaxs::hist::detail::CompactCoordinatesXYZW<vbw>::evaluate_rounded_4_sse_into(
        std::span<const CompactCoordinatesXYZW, 4> others,
        int32_t* dist_out, float* wt_out
    ) const noexcept {
        const auto& v1 = others[0]; const auto& v2 = others[1]; const auto& v3 = others[2]; const auto& v4 = others[3];
        __m128 sv = _mm_load_ps(this->data.data());
        __m128 d1 = _mm_sub_ps(sv, _mm_load_ps(v1.data.data()));
        __m128 d2 = _mm_sub_ps(sv, _mm_load_ps(v2.data.data()));
        __m128 d3 = _mm_sub_ps(sv, _mm_load_ps(v3.data.data()));
        __m128 d4 = _mm_sub_ps(sv, _mm_load_ps(v4.data.data()));
        d1 = _mm_mul_ps(d1, d1);
        d2 = _mm_mul_ps(d2, d2);
        d3 = _mm_mul_ps(d3, d3);
        d4 = _mm_mul_ps(d4, d4);

        __m128 t0 = _mm_unpacklo_ps(d1, d2);
        __m128 t1 = _mm_unpackhi_ps(d1, d2);
        __m128 t2 = _mm_unpacklo_ps(d3, d4);
        __m128 t3 = _mm_unpackhi_ps(d3, d4);
        __m128 row_x = _mm_movelh_ps(t0, t2);
        __m128 row_y = _mm_movehl_ps(t2, t0);
        __m128 row_z = _mm_movelh_ps(t1, t3);

        __m128 dist2 = _mm_add_ps(_mm_add_ps(row_x, row_y), row_z);
        __m128 dist_sqrt = _mm_sqrt_ps(dist2);
        __m128 dist_binf = _mm_mul_ps(dist_sqrt, _mm_set_ps1(get_inv_width()));
        __m128i dist_bin = _mm_cvtps_epi32(dist_binf);

        __m128 w1 = _mm_set_ps1(value.w);
        __m128 w2 = _mm_set_ps(v4.value.w, v3.value.w, v2.value.w, v1.value.w);
        __m128 weights = _mm_mul_ps(w1, w2);

        _mm_store_si128(reinterpret_cast<__m128i*>(dist_out), dist_bin);
        _mm_store_ps(wt_out, weights);
    }

    template<bool vbw>
    inline ausaxs::hist::detail::xyzw::QuadEvaluatedResultRounded ausaxs::hist::detail::CompactCoordinatesXYZW<vbw>::evaluate_rounded_4_sse(
        std::span<const CompactCoordinatesXYZW, 4> others
    ) const noexcept {
        xyzw::QuadEvaluatedResultRounded result;
        evaluate_rounded_4_sse_into(others, result.distances.data(), result.weights.data());
        return result;
    }

#endif

#if defined AUSAXS_USE_AVX2
    #include <immintrin.h>

    template<bool vbw>
    inline void ausaxs::hist::detail::CompactCoordinatesXYZW<vbw>::evaluate_8_avx_into(
        std::span<const CompactCoordinatesXYZW, 8> others,
        float* dist_out, int32_t* bin_out, float* wt_out
    ) const noexcept {
        const auto& v1 = others[0];
        const float* base = v1.data.data();
        __m256 v12 = _mm256_loadu_ps(base);
        __m256 v34 = _mm256_loadu_ps(base + 8);
        __m256 v56 = _mm256_loadu_ps(base + 16);
        __m256 v78 = _mm256_loadu_ps(base + 24);

        // extract w values (position 3 in each lane) before computing differences
        __m256 wt0 = _mm256_unpackhi_ps(v12, v34);
        __m256 wt1 = _mm256_unpackhi_ps(v56, v78);
        __m256 w_all = _mm256_shuffle_ps(wt0, wt1, _MM_SHUFFLE(3,2,3,2));
        __m256 weights = _mm256_mul_ps(_mm256_set1_ps(value.w), w_all);

        // compute squared differences (squaring before transpose gives better ILP:
        // 4 independent mul ops overlap on separate execution ports)
        __m256 svv = _mm256_broadcast_ps(reinterpret_cast<const __m128*>(this->data.data()));
        __m256 d12 = _mm256_sub_ps(svv, v12);
        __m256 d34 = _mm256_sub_ps(svv, v34);
        __m256 d56 = _mm256_sub_ps(svv, v56);
        __m256 d78 = _mm256_sub_ps(svv, v78);
        d12 = _mm256_mul_ps(d12, d12);
        d34 = _mm256_mul_ps(d34, d34);
        d56 = _mm256_mul_ps(d56, d56);
        d78 = _mm256_mul_ps(d78, d78);

        // in-lane 4x4 transpose: lo lanes → atoms {1,3,5,7}, hi lanes → atoms {2,4,6,8}
        __m256 t0 = _mm256_unpacklo_ps(d12, d34);
        __m256 t1 = _mm256_unpackhi_ps(d12, d34);
        __m256 t2 = _mm256_unpacklo_ps(d56, d78);
        __m256 t3 = _mm256_unpackhi_ps(d56, d78);
        __m256 row_x = _mm256_shuffle_ps(t0, t2, _MM_SHUFFLE(1,0,1,0));
        __m256 row_y = _mm256_shuffle_ps(t0, t2, _MM_SHUFFLE(3,2,3,2));
        __m256 row_z = _mm256_shuffle_ps(t1, t3, _MM_SHUFFLE(1,0,1,0));

        __m256 dist2 = _mm256_add_ps(_mm256_add_ps(row_x, row_y), row_z);
        __m256 dist_sqrt = _mm256_sqrt_ps(dist2);
        __m256i dist_bin = _mm256_cvtps_epi32(_mm256_mul_ps(dist_sqrt, _mm256_set1_ps(get_inv_width())));

        _mm256_store_ps(dist_out, dist_sqrt);
        _mm256_store_si256(reinterpret_cast<__m256i*>(bin_out), dist_bin);
        _mm256_store_ps(wt_out, weights);
    }

    template<bool vbw>
    inline ausaxs::hist::detail::xyzw::OctoEvaluatedResult ausaxs::hist::detail::CompactCoordinatesXYZW<vbw>::evaluate_8_avx(
        std::span<const CompactCoordinatesXYZW, 8> others
    ) const noexcept {
        xyzw::OctoEvaluatedResult result;
        evaluate_8_avx_into(others, result.distances.data(), result.distance_bins.data(), result.weights.data());
        return result;
    }

    template<bool vbw>
    inline void ausaxs::hist::detail::CompactCoordinatesXYZW<vbw>::evaluate_rounded_8_avx_into(
        std::span<const CompactCoordinatesXYZW, 8> others,
        int32_t* dist_out, float* wt_out
    ) const noexcept {
        const auto& v1 = others[0];
        const float* base = v1.data.data();
        __m256 v12 = _mm256_loadu_ps(base);
        __m256 v34 = _mm256_loadu_ps(base + 8);
        __m256 v56 = _mm256_loadu_ps(base + 16);
        __m256 v78 = _mm256_loadu_ps(base + 24);

        // extract w values (position 3 in each lane) before computing differences
        __m256 wt0 = _mm256_unpackhi_ps(v12, v34);
        __m256 wt1 = _mm256_unpackhi_ps(v56, v78);
        __m256 w_all = _mm256_shuffle_ps(wt0, wt1, _MM_SHUFFLE(3,2,3,2));
        __m256 weights = _mm256_mul_ps(_mm256_set1_ps(value.w), w_all);

        __m256 svv = _mm256_broadcast_ps(reinterpret_cast<const __m128*>(this->data.data()));
        __m256 d12 = _mm256_sub_ps(svv, v12);
        __m256 d34 = _mm256_sub_ps(svv, v34);
        __m256 d56 = _mm256_sub_ps(svv, v56);
        __m256 d78 = _mm256_sub_ps(svv, v78);
        d12 = _mm256_mul_ps(d12, d12);
        d34 = _mm256_mul_ps(d34, d34);
        d56 = _mm256_mul_ps(d56, d56);
        d78 = _mm256_mul_ps(d78, d78);

        __m256 t0 = _mm256_unpacklo_ps(d12, d34);
        __m256 t1 = _mm256_unpackhi_ps(d12, d34);
        __m256 t2 = _mm256_unpacklo_ps(d56, d78);
        __m256 t3 = _mm256_unpackhi_ps(d56, d78);
        __m256 row_x = _mm256_shuffle_ps(t0, t2, _MM_SHUFFLE(1,0,1,0));
        __m256 row_y = _mm256_shuffle_ps(t0, t2, _MM_SHUFFLE(3,2,3,2));
        __m256 row_z = _mm256_shuffle_ps(t1, t3, _MM_SHUFFLE(1,0,1,0));

        __m256 dist2 = _mm256_add_ps(_mm256_add_ps(row_x, row_y), row_z);
        __m256i dist_bin = _mm256_cvtps_epi32(_mm256_mul_ps(_mm256_sqrt_ps(dist2), _mm256_set1_ps(get_inv_width())));

        _mm256_store_si256(reinterpret_cast<__m256i*>(dist_out), dist_bin);
        _mm256_store_ps(wt_out, weights);
    }

    template<bool vbw>
    inline ausaxs::hist::detail::xyzw::OctoEvaluatedResultRounded ausaxs::hist::detail::CompactCoordinatesXYZW<vbw>::evaluate_rounded_8_avx(
        std::span<const CompactCoordinatesXYZW, 8> others
    ) const noexcept {
        xyzw::OctoEvaluatedResultRounded result;
        evaluate_rounded_8_avx_into(others, result.distances.data(), result.weights.data());
        return result;
    }
#endif

#if defined AUSAXS_USE_AVX512
    #include <immintrin.h>

    template<bool vbw>
    inline ausaxs::hist::detail::xyzw::HexaEvaluatedResult ausaxs::hist::detail::CompactCoordinatesXYZW<vbw>::evaluate_16_avx512(
        std::span<const CompactCoordinatesXYZW, 16> others
    ) const noexcept {
        const float* base = others[0].data.data();
        __m512 v03   = _mm512_loadu_ps(base);
        __m512 v47   = _mm512_loadu_ps(base + 16);
        __m512 v811  = _mm512_loadu_ps(base + 32);
        __m512 v1215 = _mm512_loadu_ps(base + 48);

        // extract weights from raw data before computing differences
        const __m512i gather_w = _mm512_setr_epi32(3, 7, 11, 15, 19, 23, 27, 31, 0, 0, 0, 0, 0, 0, 0, 0);
        __m256 w_lo = _mm512_castps512_ps256(_mm512_permutex2var_ps(v03, gather_w, v47));
        __m256 w_hi = _mm512_castps512_ps256(_mm512_permutex2var_ps(v811, gather_w, v1215));
        __m512 weights = _mm512_mul_ps(_mm512_set1_ps(value.w),
            _mm512_insertf32x8(_mm512_castps256_ps512(w_lo), w_hi, 1));

        // compute squared differences (square-first for ILP)
        __m512 svv = _mm512_broadcast_f32x4(_mm_load_ps(this->data.data()));
        __m512 d03   = _mm512_sub_ps(svv, v03);
        __m512 d47   = _mm512_sub_ps(svv, v47);
        __m512 d811  = _mm512_sub_ps(svv, v811);
        __m512 d1215 = _mm512_sub_ps(svv, v1215);
        d03   = _mm512_mul_ps(d03, d03);
        d47   = _mm512_mul_ps(d47, d47);
        d811  = _mm512_mul_ps(d811, d811);
        d1215 = _mm512_mul_ps(d1215, d1215);

        // gather x², y², z² from each pair, then combine into full 512-bit vectors
        const __m512i gather_x = _mm512_setr_epi32(0, 4, 8, 12, 16, 20, 24, 28, 0, 0, 0, 0, 0, 0, 0, 0);
        const __m512i gather_y = _mm512_setr_epi32(1, 5, 9, 13, 17, 21, 25, 29, 0, 0, 0, 0, 0, 0, 0, 0);
        const __m512i gather_z = _mm512_setr_epi32(2, 6, 10, 14, 18, 22, 26, 30, 0, 0, 0, 0, 0, 0, 0, 0);

        __m256 rx_lo = _mm512_castps512_ps256(_mm512_permutex2var_ps(d03, gather_x, d47));
        __m256 rx_hi = _mm512_castps512_ps256(_mm512_permutex2var_ps(d811, gather_x, d1215));
        __m512 row_x = _mm512_insertf32x8(_mm512_castps256_ps512(rx_lo), rx_hi, 1);

        __m256 ry_lo = _mm512_castps512_ps256(_mm512_permutex2var_ps(d03, gather_y, d47));
        __m256 ry_hi = _mm512_castps512_ps256(_mm512_permutex2var_ps(d811, gather_y, d1215));
        __m512 row_y = _mm512_insertf32x8(_mm512_castps256_ps512(ry_lo), ry_hi, 1);

        __m256 rz_lo = _mm512_castps512_ps256(_mm512_permutex2var_ps(d03, gather_z, d47));
        __m256 rz_hi = _mm512_castps512_ps256(_mm512_permutex2var_ps(d811, gather_z, d1215));
        __m512 row_z = _mm512_insertf32x8(_mm512_castps256_ps512(rz_lo), rz_hi, 1);

        __m512 dist2 = _mm512_add_ps(_mm512_add_ps(row_x, row_y), row_z);
        __m512 dist_sqrt = _mm512_sqrt_ps(dist2);
        __m512i dist_bin = _mm512_cvtps_epi32(_mm512_mul_ps(dist_sqrt, _mm512_set1_ps(get_inv_width())));

        xyzw::HexaEvaluatedResult result;
        _mm512_store_ps(reinterpret_cast<float*>(result.distances.data()), dist_sqrt);
        _mm512_store_si512(reinterpret_cast<__m512i*>(result.distance_bins.data()), dist_bin);
        _mm512_store_ps(reinterpret_cast<float*>(result.weights.data()), weights);
        return result;
    }

    template<bool vbw>
    inline ausaxs::hist::detail::xyzw::HexaEvaluatedResultRounded ausaxs::hist::detail::CompactCoordinatesXYZW<vbw>::evaluate_rounded_16_avx512(
        std::span<const CompactCoordinatesXYZW, 16> others
    ) const noexcept {
        const float* base = others[0].data.data();
        __m512 v03   = _mm512_loadu_ps(base);
        __m512 v47   = _mm512_loadu_ps(base + 16);
        __m512 v811  = _mm512_loadu_ps(base + 32);
        __m512 v1215 = _mm512_loadu_ps(base + 48);

        const __m512i gather_w = _mm512_setr_epi32(3, 7, 11, 15, 19, 23, 27, 31, 0, 0, 0, 0, 0, 0, 0, 0);
        __m256 w_lo = _mm512_castps512_ps256(_mm512_permutex2var_ps(v03, gather_w, v47));
        __m256 w_hi = _mm512_castps512_ps256(_mm512_permutex2var_ps(v811, gather_w, v1215));
        __m512 weights = _mm512_mul_ps(_mm512_set1_ps(value.w),
            _mm512_insertf32x8(_mm512_castps256_ps512(w_lo), w_hi, 1));

        __m512 svv = _mm512_broadcast_f32x4(_mm_load_ps(this->data.data()));
        __m512 d03   = _mm512_sub_ps(svv, v03);
        __m512 d47   = _mm512_sub_ps(svv, v47);
        __m512 d811  = _mm512_sub_ps(svv, v811);
        __m512 d1215 = _mm512_sub_ps(svv, v1215);
        d03   = _mm512_mul_ps(d03, d03);
        d47   = _mm512_mul_ps(d47, d47);
        d811  = _mm512_mul_ps(d811, d811);
        d1215 = _mm512_mul_ps(d1215, d1215);

        const __m512i gather_x = _mm512_setr_epi32(0, 4, 8, 12, 16, 20, 24, 28, 0, 0, 0, 0, 0, 0, 0, 0);
        const __m512i gather_y = _mm512_setr_epi32(1, 5, 9, 13, 17, 21, 25, 29, 0, 0, 0, 0, 0, 0, 0, 0);
        const __m512i gather_z = _mm512_setr_epi32(2, 6, 10, 14, 18, 22, 26, 30, 0, 0, 0, 0, 0, 0, 0, 0);

        __m256 rx_lo = _mm512_castps512_ps256(_mm512_permutex2var_ps(d03, gather_x, d47));
        __m256 rx_hi = _mm512_castps512_ps256(_mm512_permutex2var_ps(d811, gather_x, d1215));
        __m512 row_x = _mm512_insertf32x8(_mm512_castps256_ps512(rx_lo), rx_hi, 1);

        __m256 ry_lo = _mm512_castps512_ps256(_mm512_permutex2var_ps(d03, gather_y, d47));
        __m256 ry_hi = _mm512_castps512_ps256(_mm512_permutex2var_ps(d811, gather_y, d1215));
        __m512 row_y = _mm512_insertf32x8(_mm512_castps256_ps512(ry_lo), ry_hi, 1);

        __m256 rz_lo = _mm512_castps512_ps256(_mm512_permutex2var_ps(d03, gather_z, d47));
        __m256 rz_hi = _mm512_castps512_ps256(_mm512_permutex2var_ps(d811, gather_z, d1215));
        __m512 row_z = _mm512_insertf32x8(_mm512_castps256_ps512(rz_lo), rz_hi, 1);

        __m512i dist_bin = _mm512_cvtps_epi32(_mm512_mul_ps(
            _mm512_sqrt_ps(_mm512_add_ps(_mm512_add_ps(row_x, row_y), row_z)),
            _mm512_set1_ps(get_inv_width())));

        xyzw::HexaEvaluatedResultRounded result;
        _mm512_store_si512(reinterpret_cast<__m512i*>(result.distances.data()), dist_bin);
        _mm512_store_ps(reinterpret_cast<float*>(result.weights.data()), weights);
        return result;
    }
#endif