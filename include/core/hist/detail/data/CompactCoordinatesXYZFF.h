// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

/**
 * @brief This file contains custom vector instructions for efficient scattering calculations.
 * 
 * The implementation is specialized for generic systems defined by a collection of [x: float, y: float, z: float, ff_type: int32] vectors.
 * This is useful for X-ray calculations with different atomic species, where each atom may have a different form factor type.
 * Note the distinction from CompactCoordinatesXYZFF, which is specialized for systems defined by a collection of [x, y, z, w] float vectors.
 *
 * The returned form factor bin is calculated as ff_bin = ff2 + ff1 * N_ff_types, where N_ff_types is the total number of form factor types.
 */

#pragma once

#include <hist/detail/data/IntrinsicMacros.h>
#include <hist/detail/data/WidthControllers.h>
#include <hist/detail/data/IntrinsicHelpers.h>
#include <math/Vector3.h>
#include <constants/Constants.h>
#include <settings/Flags.h>
#include <form_factor/FormFactorType.h>

#include <array>
#include <cstdint>
#include <span>

namespace ausaxs::hist::detail::xyzff {
    struct EvaluatedResult {
        float distance;      // The exact distance 
        int32_t distance_bin; // The distance bin index
        int32_t ff_bin;      // The form factor bin index
    };

    struct EvaluatedResultRounded {
        int32_t distance;   // The distance bin 
        int32_t ff_bin;     // The form factor bin index
    };

    struct QuadEvaluatedResult {
        std::array<float, 4> distances;       // The raw distances (for weighted bin center calculation)
        std::array<int32_t, 4> distance_bins; // The distance bin indices (for array indexing)
        std::array<int32_t, 4> ff_bins;       // The form factor bin indices
    };

    struct QuadEvaluatedResultRounded {
        std::array<int32_t, 4> distances;   // The distance bin
        std::array<int32_t, 4> ff_bins;     // The form factor bin indices
    };

    struct alignas(32) OctoEvaluatedResult {
        std::array<float, 8> distances;       // The raw distances (for weighted bin center calculation)
        std::array<int32_t, 8> distance_bins; // The distance bin indices (for array indexing)
        std::array<int32_t, 8> ff_bins;       // The form factor bin indices
    };

    struct alignas(32) OctoEvaluatedResultRounded {
        std::array<int32_t, 8> distances;   // The distance bin
        std::array<int32_t, 8> ff_bins;     // The form factor bin indices
    };

    struct alignas(32) HexaEvaluatedResult {
        std::array<float, 16> distances;       // The raw distances
        std::array<int32_t, 16> distance_bins; // The distance bin indices
        std::array<int32_t, 16> ff_bins;       // The form factor bin indices
    };

    struct alignas(32) HexaEvaluatedResultRounded {
        std::array<int32_t, 16> distances;   // The distance bin
        std::array<int32_t, 16> ff_bins;     // The form factor bin indices
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
    template<bool variable_bin_width, bool explicit_ff = false>
    class CompactCoordinatesXYZFF : public WidthController<variable_bin_width> {
        public:
            using WidthController<variable_bin_width>::get_inv_width;
            CompactCoordinatesXYZFF() noexcept = default;
            CompactCoordinatesXYZFF(const CompactCoordinatesXYZFF& other) noexcept = default;
            CompactCoordinatesXYZFF(CompactCoordinatesXYZFF&& other) noexcept = default;
            CompactCoordinatesXYZFF& operator= (const CompactCoordinatesXYZFF& other) noexcept = default;
            CompactCoordinatesXYZFF& operator= (CompactCoordinatesXYZFF&& other) noexcept = default;

            template<numeric T>
            CompactCoordinatesXYZFF(const Vector3<T>& v, int32_t ff) noexcept : value{.pos={static_cast<float>(v.x()), static_cast<float>(v.y()), static_cast<float>(v.z())}, .ff=ff} {}
            CompactCoordinatesXYZFF(const Vector3<float>& v, int32_t ff) noexcept : value{.pos=v, .ff=ff} {}

            /**
             * @brief Calculate the @a binned distance and combined ff_bin between this and a single other CompactCoordinatesXYZFF.
             */
            xyzff::EvaluatedResultRounded evaluate_rounded(const CompactCoordinatesXYZFF& other) const noexcept;
            xyzff::EvaluatedResult evaluate(const CompactCoordinatesXYZFF& other) const noexcept;

            xyzff::QuadEvaluatedResultRounded evaluate_rounded_4(std::span<const CompactCoordinatesXYZFF, 4> others) const noexcept;
            xyzff::QuadEvaluatedResult evaluate_4(std::span<const CompactCoordinatesXYZFF, 4> others) const noexcept;

            xyzff::OctoEvaluatedResultRounded evaluate_rounded_8(std::span<const CompactCoordinatesXYZFF, 8> others) const noexcept;
            xyzff::OctoEvaluatedResult evaluate_8(std::span<const CompactCoordinatesXYZFF, 8> others) const noexcept;

            xyzff::HexaEvaluatedResultRounded evaluate_rounded_16(std::span<const CompactCoordinatesXYZFF, 16> others) const noexcept;
            xyzff::HexaEvaluatedResult evaluate_16(std::span<const CompactCoordinatesXYZFF, 16> others) const noexcept;

            union {
                struct {Vector3<float> pos; int32_t ff;} value;
                std::array<float, 4> data;
            };

        protected:
            xyzff::EvaluatedResultRounded evaluate_rounded_scalar(const CompactCoordinatesXYZFF& other) const noexcept;
            xyzff::EvaluatedResult evaluate_scalar(const CompactCoordinatesXYZFF& other) const noexcept;

            xyzff::QuadEvaluatedResultRounded evaluate_rounded_4_scalar(std::span<const CompactCoordinatesXYZFF, 4> others) const noexcept;
            xyzff::QuadEvaluatedResult evaluate_4_scalar(std::span<const CompactCoordinatesXYZFF, 4> others) const noexcept;

            #if defined AUSAXS_USE_SSE2
                xyzff::QuadEvaluatedResultRounded evaluate_rounded_4_sse(std::span<const CompactCoordinatesXYZFF, 4> others) const noexcept;
                xyzff::QuadEvaluatedResult evaluate_4_sse(std::span<const CompactCoordinatesXYZFF, 4> others) const noexcept;
                void evaluate_rounded_4_sse_into(std::span<const CompactCoordinatesXYZFF, 4> others, int32_t* dist_out, int32_t* ff_out) const noexcept;
                void evaluate_4_sse_into(std::span<const CompactCoordinatesXYZFF, 4> others, float* dist_out, int32_t* bin_out, int32_t* ff_out) const noexcept;
            #endif

            #if defined AUSAXS_USE_AVX2
                xyzff::OctoEvaluatedResultRounded evaluate_rounded_8_avx(std::span<const CompactCoordinatesXYZFF, 8> others) const noexcept;
                xyzff::OctoEvaluatedResult evaluate_8_avx(std::span<const CompactCoordinatesXYZFF, 8> others) const noexcept;
                void evaluate_rounded_8_avx_into(std::span<const CompactCoordinatesXYZFF, 8> others, int32_t* dist_out, int32_t* ff_out) const noexcept;
                void evaluate_8_avx_into(std::span<const CompactCoordinatesXYZFF, 8> others, float* dist_out, int32_t* bin_out, int32_t* ff_out) const noexcept;
            #endif

            #if defined AUSAXS_USE_AVX512
                xyzff::HexaEvaluatedResultRounded evaluate_rounded_16_avx512(std::span<const CompactCoordinatesXYZFF, 16> others) const noexcept;
                xyzff::HexaEvaluatedResult evaluate_16_avx512(std::span<const CompactCoordinatesXYZFF, 16> others) const noexcept;
            #endif
    };
    static_assert(sizeof(CompactCoordinatesXYZFF<true>) == 16,              "CompactCoordinatesXYZFF is not 16 bytes. This is required for aligning SIMD instructions.");
    static_assert(std::is_trivial_v<CompactCoordinatesXYZFF<true>>,         "CompactCoordinatesXYZFF is not trivial");
    static_assert(std::is_standard_layout_v<CompactCoordinatesXYZFF<true>>, "CompactCoordinatesXYZFF is not standard layout");
    static_assert(supports_nothrow_move_v<CompactCoordinatesXYZFF<true>>,   "CompactCoordinatesXYZFF should support nothrow move semantics.");
    static_assert(sizeof(CompactCoordinatesXYZFF<false>) == 16,             "CompactCoordinatesXYZFF is not 16 bytes. This is required for aligning SIMD instructions.");
    static_assert(std::is_trivial_v<CompactCoordinatesXYZFF<false>>,        "CompactCoordinatesXYZFF is not trivial");
    static_assert(std::is_standard_layout_v<CompactCoordinatesXYZFF<false>>,"CompactCoordinatesXYZFF is not standard layout");
    static_assert(supports_nothrow_move_v<CompactCoordinatesXYZFF<false>>,  "CompactCoordinatesXYZFF should support nothrow move semantics.");
}

//#########################################//
//############ IMPLEMENTATION #############//
//#########################################//

// implementation defined in header to support efficient inlining

namespace ausaxs::hist::detail::xyzff {
    template<bool explicit_ff>
    inline int32_t ff_bin_index(int32_t ff1, int32_t ff2) noexcept {
        if constexpr (explicit_ff) {
            return ff2 + ff1*ausaxs::form_factor::get_count_without_excluded_volume();
        } else {
            return ff2 + ff1*ausaxs::form_factor::get_count();
        }
    }
}

template<bool vbw, bool explicit_ff>
inline ausaxs::hist::detail::xyzff::EvaluatedResult ausaxs::hist::detail::CompactCoordinatesXYZFF<vbw, explicit_ff>::evaluate(const CompactCoordinatesXYZFF& other) const  noexcept{
    return evaluate_scalar(other);
}

template<bool vbw, bool explicit_ff>
inline ausaxs::hist::detail::xyzff::EvaluatedResultRounded ausaxs::hist::detail::CompactCoordinatesXYZFF<vbw, explicit_ff>::evaluate_rounded(const CompactCoordinatesXYZFF& other) const  noexcept{
    return evaluate_rounded_scalar(other);
}

template<bool vbw, bool explicit_ff>
inline ausaxs::hist::detail::xyzff::QuadEvaluatedResult ausaxs::hist::detail::CompactCoordinatesXYZFF<vbw, explicit_ff>::evaluate_4(
    std::span<const CompactCoordinatesXYZFF, 4> others
) const noexcept {
    #if defined AUSAXS_USE_SSE2
        return evaluate_4_sse(others);
    #else
        return evaluate_4_scalar(others);
    #endif
}

template<bool vbw, bool explicit_ff>
inline ausaxs::hist::detail::xyzff::QuadEvaluatedResultRounded ausaxs::hist::detail::CompactCoordinatesXYZFF<vbw, explicit_ff>::evaluate_rounded_4(
    std::span<const CompactCoordinatesXYZFF, 4> others
) const noexcept {
    #if defined AUSAXS_USE_SSE2
        return evaluate_rounded_4_sse(others);
    #else
        return evaluate_rounded_4_scalar(others);
    #endif
}

template<bool vbw, bool explicit_ff>
inline ausaxs::hist::detail::xyzff::OctoEvaluatedResult ausaxs::hist::detail::CompactCoordinatesXYZFF<vbw, explicit_ff>::evaluate_8(
    std::span<const CompactCoordinatesXYZFF, 8> others
) const noexcept {
    #if defined AUSAXS_USE_AVX2
        return evaluate_8_avx(others);
    #elif defined AUSAXS_USE_SSE2
        xyzff::OctoEvaluatedResult result;
        evaluate_4_sse_into(others.template first<4>(), result.distances.data(), result.distance_bins.data(), result.ff_bins.data());
        evaluate_4_sse_into(others.template last<4>(), result.distances.data()+4, result.distance_bins.data()+4, result.ff_bins.data()+4);
        return result;
    #else
        return xyzff::OctoEvaluatedResult(
            evaluate_4_scalar(others.template first<4>()),
            evaluate_4_scalar(others.template last<4>())
        );
    #endif
}

template<bool vbw, bool explicit_ff>
inline ausaxs::hist::detail::xyzff::OctoEvaluatedResultRounded ausaxs::hist::detail::CompactCoordinatesXYZFF<vbw, explicit_ff>::evaluate_rounded_8(
    std::span<const CompactCoordinatesXYZFF, 8> others
) const noexcept {
    #if defined AUSAXS_USE_AVX2
        return evaluate_rounded_8_avx(others);
    #elif defined AUSAXS_USE_SSE2
        xyzff::OctoEvaluatedResultRounded result;
        evaluate_rounded_4_sse_into(others.template first<4>(), result.distances.data(), result.ff_bins.data());
        evaluate_rounded_4_sse_into(others.template last<4>(), result.distances.data()+4, result.ff_bins.data()+4);
        return result;
    #else
        return xyzff::OctoEvaluatedResultRounded(
            evaluate_rounded_4_scalar(others.template first<4>()),
            evaluate_rounded_4_scalar(others.template last<4>())
        );
    #endif
}

template<bool vbw, bool explicit_ff>
inline ausaxs::hist::detail::xyzff::HexaEvaluatedResult ausaxs::hist::detail::CompactCoordinatesXYZFF<vbw, explicit_ff>::evaluate_16(
    std::span<const CompactCoordinatesXYZFF, 16> others
) const noexcept {
    #if defined AUSAXS_USE_AVX512
        return evaluate_16_avx512(others);
    #elif defined AUSAXS_USE_AVX2
        xyzff::HexaEvaluatedResult result;
        evaluate_8_avx_into(others.template first<8>(), result.distances.data(), result.distance_bins.data(), result.ff_bins.data());
        evaluate_8_avx_into(others.template last<8>(), result.distances.data()+8, result.distance_bins.data()+8, result.ff_bins.data()+8);
        return result;
    #elif defined AUSAXS_USE_SSE2
        xyzff::HexaEvaluatedResult result;
        evaluate_4_sse_into(others.template first<4>(), result.distances.data(), result.distance_bins.data(), result.ff_bins.data());
        evaluate_4_sse_into(others.template subspan<4,4>(), result.distances.data()+4, result.distance_bins.data()+4, result.ff_bins.data()+4);
        evaluate_4_sse_into(others.template subspan<8,4>(), result.distances.data()+8, result.distance_bins.data()+8, result.ff_bins.data()+8);
        evaluate_4_sse_into(others.template last<4>(), result.distances.data()+12, result.distance_bins.data()+12, result.ff_bins.data()+12);
        return result;
    #else
        return xyzff::HexaEvaluatedResult(evaluate_8(others.template first<8>()), evaluate_8(others.template last<8>()));
    #endif
}

template<bool vbw, bool explicit_ff>
inline ausaxs::hist::detail::xyzff::HexaEvaluatedResultRounded ausaxs::hist::detail::CompactCoordinatesXYZFF<vbw, explicit_ff>::evaluate_rounded_16(
    std::span<const CompactCoordinatesXYZFF, 16> others
) const noexcept {
    #if defined AUSAXS_USE_AVX512
        return evaluate_rounded_16_avx512(others);
    #elif defined AUSAXS_USE_AVX2
        xyzff::HexaEvaluatedResultRounded result;
        evaluate_rounded_8_avx_into(others.template first<8>(), result.distances.data(), result.ff_bins.data());
        evaluate_rounded_8_avx_into(others.template last<8>(), result.distances.data()+8, result.ff_bins.data()+8);
        return result;
    #elif defined AUSAXS_USE_SSE2
        xyzff::HexaEvaluatedResultRounded result;
        evaluate_rounded_4_sse_into(others.template first<4>(), result.distances.data(), result.ff_bins.data());
        evaluate_rounded_4_sse_into(others.template subspan<4,4>(), result.distances.data()+4, result.ff_bins.data()+4);
        evaluate_rounded_4_sse_into(others.template subspan<8,4>(), result.distances.data()+8, result.ff_bins.data()+8);
        evaluate_rounded_4_sse_into(others.template last<4>(), result.distances.data()+12, result.ff_bins.data()+12);
        return result;
    #else
        return xyzff::HexaEvaluatedResultRounded(evaluate_rounded_8(others.template first<8>()), evaluate_rounded_8(others.template last<8>()));
    #endif
}

template<bool vbw, bool explicit_ff>
inline ausaxs::hist::detail::xyzff::EvaluatedResult ausaxs::hist::detail::CompactCoordinatesXYZFF<vbw, explicit_ff>::evaluate_scalar(
    const CompactCoordinatesXYZFF& other
) const noexcept {
    float dist = std::sqrt(squared_dot_product(this->data.data(), other.data.data()));
    int32_t dist_bin = std::round(get_inv_width() * dist);
    int32_t ff_bin = xyzff::ff_bin_index<explicit_ff>(this->value.ff, other.value.ff);
    return xyzff::EvaluatedResult(dist, dist_bin, ff_bin);
}

template<bool vbw, bool explicit_ff>
inline ausaxs::hist::detail::xyzff::EvaluatedResultRounded ausaxs::hist::detail::CompactCoordinatesXYZFF<vbw, explicit_ff>::evaluate_rounded_scalar(
    const CompactCoordinatesXYZFF& other
) const noexcept {
    int32_t dist = std::round(get_inv_width()*std::sqrt(squared_dot_product(this->data.data(), other.data.data())));
    int32_t ff_bin = xyzff::ff_bin_index<explicit_ff>(this->value.ff, other.value.ff);
    return xyzff::EvaluatedResultRounded(dist, ff_bin);
}

template<bool vbw, bool explicit_ff>
inline ausaxs::hist::detail::xyzff::QuadEvaluatedResult ausaxs::hist::detail::CompactCoordinatesXYZFF<vbw, explicit_ff>::evaluate_4_scalar(
    std::span<const CompactCoordinatesXYZFF, 4> others
) const noexcept {
    const auto& v1 = others[0]; const auto& v2 = others[1]; const auto& v3 = others[2]; const auto& v4 = others[3];
    float dx1 = std::sqrt(squared_dot_product(this->data.data(), v1.data.data()));
    float dx2 = std::sqrt(squared_dot_product(this->data.data(), v2.data.data()));
    float dx3 = std::sqrt(squared_dot_product(this->data.data(), v3.data.data()));
    float dx4 = std::sqrt(squared_dot_product(this->data.data(), v4.data.data()));
    float inv_width = get_inv_width();
    return xyzff::QuadEvaluatedResult(
        std::array<float, 4>{dx1, dx2, dx3, dx4},
        std::array<int32_t, 4>{
            static_cast<int32_t>(std::round(inv_width * dx1)),
            static_cast<int32_t>(std::round(inv_width * dx2)),
            static_cast<int32_t>(std::round(inv_width * dx3)),
            static_cast<int32_t>(std::round(inv_width * dx4))
        },
        std::array<int32_t, 4>{
            xyzff::ff_bin_index<explicit_ff>(this->value.ff, v1.value.ff),
            xyzff::ff_bin_index<explicit_ff>(this->value.ff, v2.value.ff),
            xyzff::ff_bin_index<explicit_ff>(this->value.ff, v3.value.ff),
            xyzff::ff_bin_index<explicit_ff>(this->value.ff, v4.value.ff)
        }
    );
}

template<bool vbw, bool explicit_ff>
inline ausaxs::hist::detail::xyzff::QuadEvaluatedResultRounded ausaxs::hist::detail::CompactCoordinatesXYZFF<vbw, explicit_ff>::evaluate_rounded_4_scalar(
    std::span<const CompactCoordinatesXYZFF, 4> others
) const noexcept {
    const auto& v1 = others[0]; const auto& v2 = others[1]; const auto& v3 = others[2]; const auto& v4 = others[3];
    int32_t dx1 = std::round(get_inv_width()*std::sqrt(squared_dot_product(this->data.data(), v1.data.data())));
    int32_t dx2 = std::round(get_inv_width()*std::sqrt(squared_dot_product(this->data.data(), v2.data.data())));
    int32_t dx3 = std::round(get_inv_width()*std::sqrt(squared_dot_product(this->data.data(), v3.data.data())));
    int32_t dx4 = std::round(get_inv_width()*std::sqrt(squared_dot_product(this->data.data(), v4.data.data())));
    return xyzff::QuadEvaluatedResultRounded(
        std::array<int32_t, 4>{dx1, dx2, dx3, dx4},
        std::array<int32_t, 4>{
            xyzff::ff_bin_index<explicit_ff>(this->value.ff, v1.value.ff),
            xyzff::ff_bin_index<explicit_ff>(this->value.ff, v2.value.ff),
            xyzff::ff_bin_index<explicit_ff>(this->value.ff, v3.value.ff),
            xyzff::ff_bin_index<explicit_ff>(this->value.ff, v4.value.ff)
        }
    );
}

#if defined AUSAXS_USE_SSE2
    #include <nmmintrin.h>
    namespace ausaxs::hist::detail::xyzff {
        template<bool vbw, bool explicit_ff>
        inline static __m128i ff_bin_index(int32_t ff_self, __m128 ff_others_raw) noexcept {
            __m128 mul_fac;
            if constexpr (explicit_ff) {
                mul_fac = _mm_set_ps1(form_factor::get_count_without_excluded_volume());
            } else {
                mul_fac = _mm_set_ps1(form_factor::get_count());
            }
            __m128 ff_others = _mm_cvtepi32_ps(_mm_castps_si128(ff_others_raw));
            __m128 ff1_scaled = _mm_mul_ps(_mm_set_ps1(static_cast<float>(ff_self)), mul_fac);
            return _mm_cvtps_epi32(_mm_add_ps(ff_others, ff1_scaled));
        }
    }

    template<bool vbw, bool explicit_ff>
    inline void ausaxs::hist::detail::CompactCoordinatesXYZFF<vbw, explicit_ff>::evaluate_4_sse_into(
        std::span<const CompactCoordinatesXYZFF, 4> others,
        float* dist_out, int32_t* bin_out, int32_t* ff_out
    ) const noexcept {
        const float* p = reinterpret_cast<const float*>(others.data());
        __m128 r1 = _mm_loadu_ps(p);
        __m128 r2 = _mm_loadu_ps(p + 4);
        __m128 r3 = _mm_loadu_ps(p + 8);
        __m128 r4 = _mm_loadu_ps(p + 12);

        // extract ff values (position 3, stored as int32) before computing differences
        __m128i ff_bins = xyzff::ff_bin_index<vbw, explicit_ff>(this->value.ff,
            _mm_movehl_ps(_mm_unpackhi_ps(r3, r4), _mm_unpackhi_ps(r1, r2)));

        __m128 sv = _mm_load_ps(this->data.data());
        __m128 d1 = _mm_sub_ps(sv, r1);
        __m128 d2 = _mm_sub_ps(sv, r2);
        __m128 d3 = _mm_sub_ps(sv, r3);
        __m128 d4 = _mm_sub_ps(sv, r4);
        d1 = _mm_mul_ps(d1, d1);
        d2 = _mm_mul_ps(d2, d2);
        d3 = _mm_mul_ps(d3, d3);
        d4 = _mm_mul_ps(d4, d4);

        // 4x4 transpose: converts per-atom [dx²,dy²,dz²,dff²] to per-component rows, discarding ff
        __m128 t0 = _mm_unpacklo_ps(d1, d2);   // [dx1²,dx2²,dy1²,dy2²]
        __m128 t1 = _mm_unpackhi_ps(d1, d2);   // [dz1²,dz2²,dff1²,dff2²]
        __m128 t2 = _mm_unpacklo_ps(d3, d4);   // [dx3²,dx4²,dy3²,dy4²]
        __m128 t3 = _mm_unpackhi_ps(d3, d4);   // [dz3²,dz4²,dff3²,dff4²]
        __m128 row_x = _mm_movelh_ps(t0, t2);  // [dx1²,dx2²,dx3²,dx4²]
        __m128 row_y = _mm_movehl_ps(t2, t0);  // [dy1²,dy2²,dy3²,dy4²]
        __m128 row_z = _mm_movelh_ps(t1, t3);  // [dz1²,dz2²,dz3²,dz4²]

        __m128 dist2 = _mm_add_ps(_mm_add_ps(row_x, row_y), row_z);
        __m128 dist_sqrt = _mm_sqrt_ps(dist2);
        __m128 dist_binf = _mm_mul_ps(dist_sqrt, _mm_set_ps1(get_inv_width()));
        __m128i dist_bin = _mm_cvtps_epi32(dist_binf);

        _mm_store_ps(dist_out, dist_sqrt);
        _mm_store_si128(reinterpret_cast<__m128i*>(bin_out), dist_bin);
        _mm_store_si128(reinterpret_cast<__m128i*>(ff_out), ff_bins);
    }

    template<bool vbw, bool explicit_ff>
    inline ausaxs::hist::detail::xyzff::QuadEvaluatedResult ausaxs::hist::detail::CompactCoordinatesXYZFF<vbw, explicit_ff>::evaluate_4_sse(
        std::span<const CompactCoordinatesXYZFF, 4> others
    ) const noexcept {
        xyzff::QuadEvaluatedResult result;
        evaluate_4_sse_into(others, result.distances.data(), result.distance_bins.data(), result.ff_bins.data());
        return result;
    }

    template<bool vbw, bool explicit_ff>
    inline void ausaxs::hist::detail::CompactCoordinatesXYZFF<vbw, explicit_ff>::evaluate_rounded_4_sse_into(
        std::span<const CompactCoordinatesXYZFF, 4> others,
        int32_t* dist_out, int32_t* ff_out
    ) const noexcept {
        const float* p = reinterpret_cast<const float*>(others.data());
        __m128 r1 = _mm_loadu_ps(p);
        __m128 r2 = _mm_loadu_ps(p + 4);
        __m128 r3 = _mm_loadu_ps(p + 8);
        __m128 r4 = _mm_loadu_ps(p + 12);

        // extract ff values (position 3, stored as int32) before computing differences
        __m128i ff_bins = xyzff::ff_bin_index<vbw, explicit_ff>(this->value.ff,
            _mm_movehl_ps(_mm_unpackhi_ps(r3, r4), _mm_unpackhi_ps(r1, r2)));

        __m128 sv = _mm_load_ps(this->data.data());
        __m128 d1 = _mm_sub_ps(sv, r1);
        __m128 d2 = _mm_sub_ps(sv, r2);
        __m128 d3 = _mm_sub_ps(sv, r3);
        __m128 d4 = _mm_sub_ps(sv, r4);
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

        _mm_store_si128(reinterpret_cast<__m128i*>(dist_out), dist_bin);
        _mm_store_si128(reinterpret_cast<__m128i*>(ff_out), ff_bins);
    }

    template<bool vbw, bool explicit_ff>
    inline ausaxs::hist::detail::xyzff::QuadEvaluatedResultRounded ausaxs::hist::detail::CompactCoordinatesXYZFF<vbw, explicit_ff>::evaluate_rounded_4_sse(
        std::span<const CompactCoordinatesXYZFF, 4> others
    ) const noexcept {
        xyzff::QuadEvaluatedResultRounded result;
        evaluate_rounded_4_sse_into(others, result.distances.data(), result.ff_bins.data());
        return result;
    }

#endif

#if defined AUSAXS_USE_AVX2
    #include <immintrin.h>
    template<bool vbw, bool explicit_ff>
    inline void ausaxs::hist::detail::CompactCoordinatesXYZFF<vbw, explicit_ff>::evaluate_8_avx_into(
        std::span<const CompactCoordinatesXYZFF, 8> others,
        float* dist_out, int32_t* bin_out, int32_t* ff_out
    ) const noexcept {
        const float* p = reinterpret_cast<const float*>(others.data());
        __m256 v12 = _mm256_loadu_ps(p);
        __m256 v34 = _mm256_loadu_ps(p + 8);
        __m256 v56 = _mm256_loadu_ps(p + 16);
        __m256 v78 = _mm256_loadu_ps(p + 24);

        // extract ff values (int32 at position 3 in each lane) before computing differences
        __m256 fft0 = _mm256_unpackhi_ps(v12, v34);
        __m256 fft1 = _mm256_unpackhi_ps(v56, v78);
        __m256 ff_raw = _mm256_shuffle_ps(fft0, fft1, _MM_SHUFFLE(3,2,3,2));
        __m256 ff_float = _mm256_cvtepi32_ps(_mm256_castps_si256(ff_raw));
        __m256 mul_fac;
        if constexpr (explicit_ff) {
            mul_fac = _mm256_set1_ps(form_factor::get_count_without_excluded_volume());
        } else {
            mul_fac = _mm256_set1_ps(form_factor::get_count());
        }
        __m256i ff_bins = _mm256_cvtps_epi32(_mm256_add_ps(ff_float, _mm256_mul_ps(_mm256_set1_ps(this->value.ff), mul_fac)));

        // compute differences and square
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
        _mm256_store_si256(reinterpret_cast<__m256i*>(ff_out), ff_bins);
    }

    template<bool vbw, bool explicit_ff>
    inline ausaxs::hist::detail::xyzff::OctoEvaluatedResult ausaxs::hist::detail::CompactCoordinatesXYZFF<vbw, explicit_ff>::evaluate_8_avx(
        std::span<const CompactCoordinatesXYZFF, 8> others
    ) const noexcept {
        xyzff::OctoEvaluatedResult result;
        evaluate_8_avx_into(others, result.distances.data(), result.distance_bins.data(), result.ff_bins.data());
        return result;
    }

    template<bool vbw, bool explicit_ff>
    inline void ausaxs::hist::detail::CompactCoordinatesXYZFF<vbw, explicit_ff>::evaluate_rounded_8_avx_into(
        std::span<const CompactCoordinatesXYZFF, 8> others,
        int32_t* dist_out, int32_t* ff_out
    ) const noexcept {
        const float* p = reinterpret_cast<const float*>(others.data());
        __m256 v12 = _mm256_loadu_ps(p);
        __m256 v34 = _mm256_loadu_ps(p + 8);
        __m256 v56 = _mm256_loadu_ps(p + 16);
        __m256 v78 = _mm256_loadu_ps(p + 24);

        __m256 fft0 = _mm256_unpackhi_ps(v12, v34);
        __m256 fft1 = _mm256_unpackhi_ps(v56, v78);
        __m256 ff_raw = _mm256_shuffle_ps(fft0, fft1, _MM_SHUFFLE(3,2,3,2));
        __m256 ff_float = _mm256_cvtepi32_ps(_mm256_castps_si256(ff_raw));
        __m256 mul_fac;
        if constexpr (explicit_ff) {
            mul_fac = _mm256_set1_ps(form_factor::get_count_without_excluded_volume());
        } else {
            mul_fac = _mm256_set1_ps(form_factor::get_count());
        }
        __m256i ff_bins = _mm256_cvtps_epi32(_mm256_add_ps(ff_float, _mm256_mul_ps(_mm256_set1_ps(this->value.ff), mul_fac)));

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

        __m256i dist_bin = _mm256_cvtps_epi32(_mm256_mul_ps(_mm256_sqrt_ps(_mm256_add_ps(_mm256_add_ps(row_x, row_y), row_z)), _mm256_set1_ps(get_inv_width())));

        _mm256_store_si256(reinterpret_cast<__m256i*>(dist_out), dist_bin);
        _mm256_store_si256(reinterpret_cast<__m256i*>(ff_out), ff_bins);
    }

    template<bool vbw, bool explicit_ff>
    inline ausaxs::hist::detail::xyzff::OctoEvaluatedResultRounded ausaxs::hist::detail::CompactCoordinatesXYZFF<vbw, explicit_ff>::evaluate_rounded_8_avx(
        std::span<const CompactCoordinatesXYZFF, 8> others
    ) const noexcept {
        xyzff::OctoEvaluatedResultRounded result;
        evaluate_rounded_8_avx_into(others, result.distances.data(), result.ff_bins.data());
        return result;
    }
#endif

#if defined AUSAXS_USE_AVX512
    #include <immintrin.h>

    template<bool vbw, bool explicit_ff>
    inline ausaxs::hist::detail::xyzff::HexaEvaluatedResult ausaxs::hist::detail::CompactCoordinatesXYZFF<vbw, explicit_ff>::evaluate_16_avx512(
        std::span<const CompactCoordinatesXYZFF, 16> others
    ) const noexcept {
        const float* p = reinterpret_cast<float*>(others.data());
        __m512 v03   = _mm512_loadu_ps(p);
        __m512 v47   = _mm512_loadu_ps(p + 16);
        __m512 v811  = _mm512_loadu_ps(p + 32);
        __m512 v1215 = _mm512_loadu_ps(p + 48);

        // extract ff values from raw data before computing differences
        const __m512i gather_ff = _mm512_setr_epi32(3, 7, 11, 15, 19, 23, 27, 31, 0, 0, 0, 0, 0, 0, 0, 0);
        __m256 ff_lo_raw = _mm512_castps512_ps256(_mm512_permutex2var_ps(v03, gather_ff, v47));
        __m256 ff_hi_raw = _mm512_castps512_ps256(_mm512_permutex2var_ps(v811, gather_ff, v1215));
        __m512 ff_all_float = _mm512_cvtepi32_ps(_mm512_castps_si512(
            _mm512_insertf32x8(_mm512_castps256_ps512(ff_lo_raw), ff_hi_raw, 1)));
        __m512 mul_fac;
        if constexpr (explicit_ff) {
            mul_fac = _mm512_set1_ps(form_factor::get_count_without_excluded_volume());
        } else {
            mul_fac = _mm512_set1_ps(form_factor::get_count());
        }
        __m512i ff_bins = _mm512_cvtps_epi32(_mm512_add_ps(ff_all_float, _mm512_mul_ps(_mm512_set1_ps(this->value.ff), mul_fac)));

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

        xyzff::HexaEvaluatedResult result;
        _mm512_store_ps(reinterpret_cast<float*>(result.distances.data()), dist_sqrt);
        _mm512_store_si512(reinterpret_cast<__m512i*>(result.distance_bins.data()), dist_bin);
        _mm512_store_si512(reinterpret_cast<__m512i*>(result.ff_bins.data()), ff_bins);
        return result;
    }

    template<bool vbw, bool explicit_ff>
    inline ausaxs::hist::detail::xyzff::HexaEvaluatedResultRounded ausaxs::hist::detail::CompactCoordinatesXYZFF<vbw, explicit_ff>::evaluate_rounded_16_avx512(
        std::span<const CompactCoordinatesXYZFF, 16> others
    ) const noexcept {
        const float* p = reinterpret_cast<float*>(others.data());
        __m512 v03   = _mm512_loadu_ps(p);
        __m512 v47   = _mm512_loadu_ps(p + 16);
        __m512 v811  = _mm512_loadu_ps(p + 32);
        __m512 v1215 = _mm512_loadu_ps(p + 48);

        const __m512i gather_ff = _mm512_setr_epi32(3, 7, 11, 15, 19, 23, 27, 31, 0, 0, 0, 0, 0, 0, 0, 0);
        __m256 ff_lo_raw = _mm512_castps512_ps256(_mm512_permutex2var_ps(v03, gather_ff, v47));
        __m256 ff_hi_raw = _mm512_castps512_ps256(_mm512_permutex2var_ps(v811, gather_ff, v1215));
        __m512 ff_all_float = _mm512_cvtepi32_ps(_mm512_castps_si512(
            _mm512_insertf32x8(_mm512_castps256_ps512(ff_lo_raw), ff_hi_raw, 1)));
        __m512 mul_fac;
        if constexpr (explicit_ff) {
            mul_fac = _mm512_set1_ps(form_factor::get_count_without_excluded_volume());
        } else {
            mul_fac = _mm512_set1_ps(form_factor::get_count());
        }
        __m512i ff_bins = _mm512_cvtps_epi32(_mm512_add_ps(ff_all_float, _mm512_mul_ps(_mm512_set1_ps(this->value.ff), mul_fac)));

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

        xyzff::HexaEvaluatedResultRounded result;
        _mm512_store_si512(reinterpret_cast<__m512i*>(result.distances.data()), dist_bin);
        _mm512_store_si512(reinterpret_cast<__m512i*>(result.ff_bins.data()), ff_bins);
        return result;
    }
#endif