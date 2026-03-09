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

namespace ausaxs::hist::detail::xyzff {
    /**
     * @brief Simple structure for storing the results of a distance and form factor bin calculation.
     */
    struct EvaluatedResult {
        EvaluatedResult() noexcept = default;
        EvaluatedResult(float distance, int32_t distance_bin, int32_t ff_bin) noexcept : distance(distance), distance_bin(distance_bin), ff_bin(ff_bin) {}
        float distance;      // The exact distance 
        int32_t distance_bin; // The distance bin index
        int32_t ff_bin;      // The form factor bin index
    };

    struct EvaluatedResultRounded {
        EvaluatedResultRounded() noexcept = default;
        EvaluatedResultRounded(int32_t distance_bin, float ff_bin) noexcept : distance(distance_bin), ff_bin(ff_bin) {}
        int32_t distance;   // The distance bin 
        int32_t ff_bin;     // The form factor bin index
    };

    /**
     * @brief Simple structure for storing the results of four distance and form factor bin calculations.
     *        This is necessary to restructure the memory layout for more efficient 128-bit SIMD instructions.
     */
    struct QuadEvaluatedResult {
        QuadEvaluatedResult() noexcept = default;
        QuadEvaluatedResult(const EvaluatedResult& v1, const EvaluatedResult& v2, const EvaluatedResult& v3, const EvaluatedResult& v4) noexcept
            : distances{v1.distance, v2.distance, v3.distance, v4.distance}, 
              distance_bins{v1.distance_bin, v2.distance_bin, v3.distance_bin, v4.distance_bin},
              ff_bins{v1.ff_bin, v2.ff_bin, v3.ff_bin, v4.ff_bin}
        {}
        QuadEvaluatedResult(const std::array<float, 4>& distances, const std::array<int32_t, 4>& distance_bins, const std::array<int32_t, 4>& ff_bins) noexcept 
            : distances(distances), distance_bins(distance_bins), ff_bins(ff_bins) 
        {}

        std::array<float, 4> distances;       // The raw distances (for weighted bin center calculation)
        std::array<int32_t, 4> distance_bins; // The distance bin indices (for array indexing)
        std::array<int32_t, 4> ff_bins;       // The form factor bin indices
    };

    /**
     * @brief Simple structure for storing the results of four distance and form factor bin calculations.
     *        This is necessary to restructure the memory layout for more efficient 128-bit SIMD instructions.
     */
    struct QuadEvaluatedResultRounded {
        QuadEvaluatedResultRounded() noexcept = default;
        QuadEvaluatedResultRounded(const EvaluatedResultRounded& v1, const EvaluatedResultRounded& v2, const EvaluatedResultRounded& v3, const EvaluatedResultRounded& v4) noexcept 
            : distances{v1.distance, v2.distance, v3.distance, v4.distance}, 
                ff_bins{v1.ff_bin, v2.ff_bin, v3.ff_bin, v4.ff_bin} 
        {}
        QuadEvaluatedResultRounded(const std::array<int32_t, 4>& distances, const std::array<int32_t, 4>& ff_bins) noexcept 
            : distances(distances), ff_bins(ff_bins) 
        {}

        std::array<int32_t, 4> distances;   // The distance bin
        std::array<int32_t, 4> ff_bins;     // The form factor bin indices
    };

    /**
     * @brief Simple structure for storing the results of eight distance and form factor bin calculations.
     *        This is necessary to restructure the memory layout for more efficient 256-bit SIMD instructions.
     */
    struct alignas(32) OctoEvaluatedResult {
        OctoEvaluatedResult() noexcept = default;
        OctoEvaluatedResult(
            const EvaluatedResult& v1, const EvaluatedResult& v2, const EvaluatedResult& v3, const EvaluatedResult& v4, 
            const EvaluatedResult& v5, const EvaluatedResult& v6, const EvaluatedResult& v7, const EvaluatedResult& v8) noexcept 
            : distances{v1.distance, v2.distance, v3.distance, v4.distance, v5.distance, v6.distance, v7.distance, v8.distance},
              distance_bins{v1.distance_bin, v2.distance_bin, v3.distance_bin, v4.distance_bin, v5.distance_bin, v6.distance_bin, v7.distance_bin, v8.distance_bin},
              ff_bins{v1.ff_bin, v2.ff_bin, v3.ff_bin, v4.ff_bin, v5.ff_bin, v6.ff_bin, v7.ff_bin, v8.ff_bin}
        {}
        OctoEvaluatedResult(const std::array<float, 8>& distances, const std::array<int32_t, 8>& distance_bins, const std::array<int32_t, 8>& ff_bins) noexcept 
            : distances(distances), distance_bins(distance_bins), ff_bins(ff_bins) 
        {}

        std::array<float, 8> distances;       // The raw distances (for weighted bin center calculation)
        std::array<int32_t, 8> distance_bins; // The distance bin indices (for array indexing)
        std::array<int32_t, 8> ff_bins;       // The form factor bin indices
    };

    /**
     * @brief Simple structure for storing the results of eight distance and form factor bin calculations.
     *        This is necessary to restructure the memory layout for more efficient 256-bit SIMD instructions.
     */
    struct alignas(32) OctoEvaluatedResultRounded {
        OctoEvaluatedResultRounded() noexcept = default;
        OctoEvaluatedResultRounded(
            const EvaluatedResultRounded& v1, const EvaluatedResultRounded& v2, const EvaluatedResultRounded& v3, const EvaluatedResultRounded& v4, 
            const EvaluatedResultRounded& v5, const EvaluatedResultRounded& v6, const EvaluatedResultRounded& v7, const EvaluatedResultRounded& v8) noexcept 
            : distances{v1.distance, v2.distance, v3.distance, v4.distance, v5.distance, v6.distance, v7.distance, v8.distance}, 
              ff_bins{v1.ff_bin, v2.ff_bin, v3.ff_bin, v4.ff_bin, v5.ff_bin, v6.ff_bin, v7.ff_bin, v8.ff_bin} 
        {}
        OctoEvaluatedResultRounded(const std::array<int32_t, 8>& distances, const std::array<int32_t, 8>& ff_bins) noexcept 
            : distances(distances), ff_bins(ff_bins) 
        {}

        std::array<int32_t, 8> distances;   // The distance bin
        std::array<int32_t, 8> ff_bins;     // The form factor bin indices
    };

    // assert that it is safe to perform memcpy and reinterpret_cast on these structures
    static_assert(sizeof(EvaluatedResult)            == 12, "hist::detail::EvaluatedResult is not 12 bytes long");
    static_assert(sizeof(EvaluatedResultRounded)     == 8,  "hist::detail::EvaluatedResultRounded is not 8 bytes long");
    static_assert(sizeof(QuadEvaluatedResult)        == 48, "hist::detail::QuadEvaluatedResult is not 48 bytes long");
    static_assert(sizeof(QuadEvaluatedResultRounded) == 32, "hist::detail::QuadEvaluatedResultRounded is not 32 bytes long");
    static_assert(sizeof(OctoEvaluatedResult)        == 96, "hist::detail::OctoEvaluatedResult is not 96 bytes long");
    static_assert(sizeof(OctoEvaluatedResultRounded) == 64, "hist::detail::OctoEvaluatedResultRounded is not 64 bytes long");

    // ensure our structures are trivially copyable
    static_assert(std::is_trivial_v<EvaluatedResult>,            "hist::detail::EvaluatedResult is not trivial");
    static_assert(std::is_trivial_v<EvaluatedResultRounded>,     "hist::detail::EvaluatedResultRounded is not trivial");
    static_assert(std::is_trivial_v<QuadEvaluatedResult>,        "hist::detail::QuadEvaluatedResult is not trivial");
    static_assert(std::is_trivial_v<QuadEvaluatedResultRounded>, "hist::detail::QuadEvaluatedResultRounded is not trivial");
    static_assert(std::is_trivial_v<OctoEvaluatedResult>,        "hist::detail::OctoEvaluatedResult is not trivial");
    static_assert(std::is_trivial_v<OctoEvaluatedResultRounded>, "hist::detail::OctoEvaluatedResultRounded is not trivial");

    // check that the structures have a standard memory layout. this is required for the reinterpret_casts.
    static_assert(std::is_standard_layout_v<EvaluatedResult>,            "hist::detail::EvaluatedResult is not trivial");
    static_assert(std::is_standard_layout_v<EvaluatedResultRounded>,     "hist::detail::EvaluatedResultRounded is not trivial");
    static_assert(std::is_standard_layout_v<QuadEvaluatedResult>,        "hist::detail::QuadEvaluatedResult is not trivial");
    static_assert(std::is_standard_layout_v<QuadEvaluatedResultRounded>, "hist::detail::QuadEvaluatedResultRounded is not trivial");
    static_assert(std::is_standard_layout_v<OctoEvaluatedResult>,        "hist::detail::OctoEvaluatedResult is not trivial");
    static_assert(std::is_standard_layout_v<OctoEvaluatedResultRounded>, "hist::detail::OctoEvaluatedResultRounded is not trivial");
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

            /**
             * @brief Calculate the distance and combined ff_bin between this and a single other CompactCoordinatesXYZFF.
             */
            xyzff::EvaluatedResult evaluate(const CompactCoordinatesXYZFF& other) const noexcept;

            /**
             * @brief Calculate the @a binned distance and combined weight between this and four other CompactCoordinatesXYZFF.
             *        This leverages more efficient SIMD instructions by using a 128-bit registers (SSE).
             */
            xyzff::QuadEvaluatedResultRounded evaluate_rounded(const CompactCoordinatesXYZFF& v1, const CompactCoordinatesXYZFF& v2, const CompactCoordinatesXYZFF& v3, const CompactCoordinatesXYZFF& v4) const noexcept;

            /**
             * @brief Calculate the distance and combined weight between this and four other CompactCoordinatesXYZFF.
             *        This leverages more efficient SIMD instructions by using a 128-bit registers (SSE).
             */
            xyzff::QuadEvaluatedResult evaluate(const CompactCoordinatesXYZFF& v1, const CompactCoordinatesXYZFF& v2, const CompactCoordinatesXYZFF& v3, const CompactCoordinatesXYZFF& v4) const noexcept;

            /**
             * @brief Calculate the @a binned distance and combined weight between this and four other CompactCoordinatesXYZFF.
             *        This leverages more efficient SIMD instructions by using either two 128-bit registers (SSE) or one 256-bit register (AVX).
             */
            xyzff::OctoEvaluatedResultRounded evaluate_rounded(
                const CompactCoordinatesXYZFF& v1, const CompactCoordinatesXYZFF& v2, const CompactCoordinatesXYZFF& v3, const CompactCoordinatesXYZFF& v4,
                const CompactCoordinatesXYZFF& v5, const CompactCoordinatesXYZFF& v6, const CompactCoordinatesXYZFF& v7, const CompactCoordinatesXYZFF& v8
            ) const noexcept;

            /**
             * @brief Calculate the distance and combined weight between this and four other CompactCoordinatesXYZFF.
             *        This leverages more efficient SIMD instructions by using either two 128-bit registers (SSE) or one 256-bit register (AVX).
             */
            xyzff::OctoEvaluatedResult evaluate(
                const CompactCoordinatesXYZFF& v1, const CompactCoordinatesXYZFF& v2, const CompactCoordinatesXYZFF& v3, const CompactCoordinatesXYZFF& v4,
                const CompactCoordinatesXYZFF& v5, const CompactCoordinatesXYZFF& v6, const CompactCoordinatesXYZFF& v7, const CompactCoordinatesXYZFF& v8
            ) const noexcept;

            union {
                struct {Vector3<float> pos; int32_t ff;} value;
                std::array<float, 4> data;
            };

        protected:
            xyzff::EvaluatedResultRounded evaluate_rounded_scalar(const CompactCoordinatesXYZFF& other) const noexcept;
            xyzff::EvaluatedResult evaluate_scalar(const CompactCoordinatesXYZFF& other) const noexcept;

            xyzff::QuadEvaluatedResultRounded evaluate_rounded_scalar(const CompactCoordinatesXYZFF& v1, const CompactCoordinatesXYZFF& v2, const CompactCoordinatesXYZFF& v3, const CompactCoordinatesXYZFF& v4) const noexcept;
            xyzff::QuadEvaluatedResult evaluate_scalar(const CompactCoordinatesXYZFF& v1, const CompactCoordinatesXYZFF& v2, const CompactCoordinatesXYZFF& v3, const CompactCoordinatesXYZFF& v4) const noexcept;

            xyzff::OctoEvaluatedResultRounded evaluate_rounded_scalar(
                const CompactCoordinatesXYZFF& v1, const CompactCoordinatesXYZFF& v2, const CompactCoordinatesXYZFF& v3, const CompactCoordinatesXYZFF& v4,
                const CompactCoordinatesXYZFF& v5, const CompactCoordinatesXYZFF& v6, const CompactCoordinatesXYZFF& v7, const CompactCoordinatesXYZFF& v8
            ) const noexcept;
            xyzff::OctoEvaluatedResult evaluate_scalar(
                const CompactCoordinatesXYZFF& v1, const CompactCoordinatesXYZFF& v2, const CompactCoordinatesXYZFF& v3, const CompactCoordinatesXYZFF& v4,
                const CompactCoordinatesXYZFF& v5, const CompactCoordinatesXYZFF& v6, const CompactCoordinatesXYZFF& v7, const CompactCoordinatesXYZFF& v8
            ) const noexcept;

            #if defined AUSAXS_USE_SSE2
                xyzff::EvaluatedResultRounded evaluate_rounded_sse(const CompactCoordinatesXYZFF& other) const noexcept;
                xyzff::EvaluatedResult evaluate_sse(const CompactCoordinatesXYZFF& other) const noexcept;

                xyzff::QuadEvaluatedResultRounded evaluate_rounded_sse(const CompactCoordinatesXYZFF& v1, const CompactCoordinatesXYZFF& v2, const CompactCoordinatesXYZFF& v3, const CompactCoordinatesXYZFF& v4) const noexcept;
                xyzff::QuadEvaluatedResult evaluate_sse(const CompactCoordinatesXYZFF& v1, const CompactCoordinatesXYZFF& v2, const CompactCoordinatesXYZFF& v3, const CompactCoordinatesXYZFF& v4) const noexcept;

                xyzff::OctoEvaluatedResultRounded evaluate_rounded_sse(
                    const CompactCoordinatesXYZFF& v1, const CompactCoordinatesXYZFF& v2, const CompactCoordinatesXYZFF& v3, const CompactCoordinatesXYZFF& v4,
                    const CompactCoordinatesXYZFF& v5, const CompactCoordinatesXYZFF& v6, const CompactCoordinatesXYZFF& v7, const CompactCoordinatesXYZFF& v8
                ) const noexcept;
                xyzff::OctoEvaluatedResult evaluate_sse(
                    const CompactCoordinatesXYZFF& v1, const CompactCoordinatesXYZFF& v2, const CompactCoordinatesXYZFF& v3, const CompactCoordinatesXYZFF& v4,
                    const CompactCoordinatesXYZFF& v5, const CompactCoordinatesXYZFF& v6, const CompactCoordinatesXYZFF& v7, const CompactCoordinatesXYZFF& v8
                ) const noexcept;
            #endif

            #if defined AUSAXS_USE_AVX2
                xyzff::EvaluatedResultRounded evaluate_rounded_avx(const CompactCoordinatesXYZFF& other) const noexcept;
                xyzff::EvaluatedResult evaluate_avx(const CompactCoordinatesXYZFF& other) const noexcept;

                xyzff::QuadEvaluatedResultRounded evaluate_rounded_avx(const CompactCoordinatesXYZFF& v1, const CompactCoordinatesXYZFF& v2, const CompactCoordinatesXYZFF& v3, const CompactCoordinatesXYZFF& v4) const noexcept;
                xyzff::QuadEvaluatedResult evaluate_avx(const CompactCoordinatesXYZFF& v1, const CompactCoordinatesXYZFF& v2, const CompactCoordinatesXYZFF& v3, const CompactCoordinatesXYZFF& v4) const noexcept;

                xyzff::OctoEvaluatedResultRounded evaluate_rounded_avx(
                    const CompactCoordinatesXYZFF& v1, const CompactCoordinatesXYZFF& v2, const CompactCoordinatesXYZFF& v3, const CompactCoordinatesXYZFF& v4,
                    const CompactCoordinatesXYZFF& v5, const CompactCoordinatesXYZFF& v6, const CompactCoordinatesXYZFF& v7, const CompactCoordinatesXYZFF& v8
                ) const noexcept;
                xyzff::OctoEvaluatedResult evaluate_avx(
                    const CompactCoordinatesXYZFF& v1, const CompactCoordinatesXYZFF& v2, const CompactCoordinatesXYZFF& v3, const CompactCoordinatesXYZFF& v4,
                    const CompactCoordinatesXYZFF& v5, const CompactCoordinatesXYZFF& v6, const CompactCoordinatesXYZFF& v7, const CompactCoordinatesXYZFF& v8
                ) const noexcept;
            #endif

            #if defined AUSAXS_USE_AVX512
                xyzff::OctoEvaluatedResultRounded evaluate_rounded_avx512(
                    const CompactCoordinatesXYZFF& v1, const CompactCoordinatesXYZFF& v2, const CompactCoordinatesXYZFF& v3, const CompactCoordinatesXYZFF& v4,
                    const CompactCoordinatesXYZFF& v5, const CompactCoordinatesXYZFF& v6, const CompactCoordinatesXYZFF& v7, const CompactCoordinatesXYZFF& v8
                ) const noexcept;
                xyzff::OctoEvaluatedResult evaluate_avx512(
                    const CompactCoordinatesXYZFF& v1, const CompactCoordinatesXYZFF& v2, const CompactCoordinatesXYZFF& v3, const CompactCoordinatesXYZFF& v4,
                    const CompactCoordinatesXYZFF& v5, const CompactCoordinatesXYZFF& v6, const CompactCoordinatesXYZFF& v7, const CompactCoordinatesXYZFF& v8
                ) const noexcept;
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
    #if defined AUSAXS_USE_SSE2
        return evaluate_sse(other);
    #else
        return evaluate_scalar(other);
    #endif
}

template<bool vbw, bool explicit_ff>
inline ausaxs::hist::detail::xyzff::EvaluatedResultRounded ausaxs::hist::detail::CompactCoordinatesXYZFF<vbw, explicit_ff>::evaluate_rounded(const CompactCoordinatesXYZFF& other) const  noexcept{
    #if defined AUSAXS_USE_SSE2
        return evaluate_rounded_sse(other);
    #else
        return evaluate_rounded_scalar(other);
    #endif
}

template<bool vbw, bool explicit_ff>
inline ausaxs::hist::detail::xyzff::QuadEvaluatedResult ausaxs::hist::detail::CompactCoordinatesXYZFF<vbw, explicit_ff>::evaluate(
    const CompactCoordinatesXYZFF& v1, const CompactCoordinatesXYZFF& v2, const CompactCoordinatesXYZFF& v3, const CompactCoordinatesXYZFF& v4
) const noexcept {
    #if defined AUSAXS_USE_AVX2
        return evaluate_avx(v1, v2, v3, v4);
    #elif defined AUSAXS_USE_SSE2
        return evaluate_sse(v1, v2, v3, v4);
    #else
        return evaluate_scalar(v1, v2, v3, v4);
    #endif
}

template<bool vbw, bool explicit_ff>
inline ausaxs::hist::detail::xyzff::QuadEvaluatedResultRounded ausaxs::hist::detail::CompactCoordinatesXYZFF<vbw, explicit_ff>::evaluate_rounded(
    const CompactCoordinatesXYZFF& v1, const CompactCoordinatesXYZFF& v2, const CompactCoordinatesXYZFF& v3, const CompactCoordinatesXYZFF& v4
) const noexcept {
    #if defined AUSAXS_USE_AVX2
        return evaluate_rounded_avx(v1, v2, v3, v4);
    #elif defined AUSAXS_USE_SSE2
        return evaluate_rounded_sse(v1, v2, v3, v4);
    #else
        return evaluate_rounded_scalar(v1, v2, v3, v4);
    #endif
}

template<bool vbw, bool explicit_ff>
inline ausaxs::hist::detail::xyzff::OctoEvaluatedResult ausaxs::hist::detail::CompactCoordinatesXYZFF<vbw, explicit_ff>::evaluate(
    const CompactCoordinatesXYZFF& v1, const CompactCoordinatesXYZFF& v2, const CompactCoordinatesXYZFF& v3, const CompactCoordinatesXYZFF& v4, 
    const CompactCoordinatesXYZFF& v5, const CompactCoordinatesXYZFF& v6, const CompactCoordinatesXYZFF& v7, const CompactCoordinatesXYZFF& v8
) const noexcept {
    #if defined AUSAXS_USE_AVX512
        return evaluate_avx512(v1, v2, v3, v4, v5, v6, v7, v8);
    #elif defined AUSAXS_USE_AVX2
        return evaluate_avx(v1, v2, v3, v4, v5, v6, v7, v8);
    #elif defined AUSAXS_USE_SSE2
        return evaluate_sse(v1, v2, v3, v4, v5, v6, v7, v8);
    #else
        return evaluate_scalar(v1, v2, v3, v4, v5, v6, v7, v8);
    #endif
}

template<bool vbw, bool explicit_ff>
inline ausaxs::hist::detail::xyzff::OctoEvaluatedResultRounded ausaxs::hist::detail::CompactCoordinatesXYZFF<vbw, explicit_ff>::evaluate_rounded(
    const CompactCoordinatesXYZFF& v1, const CompactCoordinatesXYZFF& v2, const CompactCoordinatesXYZFF& v3, const CompactCoordinatesXYZFF& v4, 
    const CompactCoordinatesXYZFF& v5, const CompactCoordinatesXYZFF& v6, const CompactCoordinatesXYZFF& v7, const CompactCoordinatesXYZFF& v8
) const noexcept {
    #if defined AUSAXS_USE_AVX512
        return evaluate_rounded_avx512(v1, v2, v3, v4, v5, v6, v7, v8);
    #elif defined AUSAXS_USE_AVX2
        return evaluate_rounded_avx(v1, v2, v3, v4, v5, v6, v7, v8);
    #elif defined AUSAXS_USE_SSE2
        return evaluate_rounded_sse(v1, v2, v3, v4, v5, v6, v7, v8);
    #else
        return evaluate_rounded_scalar(v1, v2, v3, v4, v5, v6, v7, v8);
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
inline ausaxs::hist::detail::xyzff::QuadEvaluatedResult ausaxs::hist::detail::CompactCoordinatesXYZFF<vbw, explicit_ff>::evaluate_scalar(
    const CompactCoordinatesXYZFF& v1, const CompactCoordinatesXYZFF& v2, const CompactCoordinatesXYZFF& v3, const CompactCoordinatesXYZFF& v4
) const noexcept {
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
inline ausaxs::hist::detail::xyzff::QuadEvaluatedResultRounded ausaxs::hist::detail::CompactCoordinatesXYZFF<vbw, explicit_ff>::evaluate_rounded_scalar(
    const CompactCoordinatesXYZFF& v1, const CompactCoordinatesXYZFF& v2, const CompactCoordinatesXYZFF& v3, const CompactCoordinatesXYZFF& v4
) const noexcept {
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

template<bool vbw, bool explicit_ff>
inline ausaxs::hist::detail::xyzff::OctoEvaluatedResult ausaxs::hist::detail::CompactCoordinatesXYZFF<vbw, explicit_ff>::evaluate_scalar(
    const CompactCoordinatesXYZFF& v1, const CompactCoordinatesXYZFF& v2, const CompactCoordinatesXYZFF& v3, const CompactCoordinatesXYZFF& v4, 
    const CompactCoordinatesXYZFF& v5, const CompactCoordinatesXYZFF& v6, const CompactCoordinatesXYZFF& v7, const CompactCoordinatesXYZFF& v8
) const noexcept {
    float dx1 = std::sqrt(squared_dot_product(this->data.data(), v1.data.data()));
    float dx2 = std::sqrt(squared_dot_product(this->data.data(), v2.data.data()));
    float dx3 = std::sqrt(squared_dot_product(this->data.data(), v3.data.data()));
    float dx4 = std::sqrt(squared_dot_product(this->data.data(), v4.data.data()));
    float dx5 = std::sqrt(squared_dot_product(this->data.data(), v5.data.data()));
    float dx6 = std::sqrt(squared_dot_product(this->data.data(), v6.data.data()));
    float dx7 = std::sqrt(squared_dot_product(this->data.data(), v7.data.data()));
    float dx8 = std::sqrt(squared_dot_product(this->data.data(), v8.data.data()));
    float inv_width = get_inv_width();
    return xyzff::OctoEvaluatedResult(
        std::array<float, 8>{dx1, dx2, dx3, dx4, dx5, dx6, dx7, dx8},
        std::array<int32_t, 8>{
            static_cast<int32_t>(std::round(inv_width * dx1)),
            static_cast<int32_t>(std::round(inv_width * dx2)),
            static_cast<int32_t>(std::round(inv_width * dx3)),
            static_cast<int32_t>(std::round(inv_width * dx4)),
            static_cast<int32_t>(std::round(inv_width * dx5)),
            static_cast<int32_t>(std::round(inv_width * dx6)),
            static_cast<int32_t>(std::round(inv_width * dx7)),
            static_cast<int32_t>(std::round(inv_width * dx8))
        },
        std::array<int32_t, 8>{
            xyzff::ff_bin_index<explicit_ff>(this->value.ff, v1.value.ff),
            xyzff::ff_bin_index<explicit_ff>(this->value.ff, v2.value.ff),
            xyzff::ff_bin_index<explicit_ff>(this->value.ff, v3.value.ff),
            xyzff::ff_bin_index<explicit_ff>(this->value.ff, v4.value.ff),
            xyzff::ff_bin_index<explicit_ff>(this->value.ff, v5.value.ff),
            xyzff::ff_bin_index<explicit_ff>(this->value.ff, v6.value.ff),
            xyzff::ff_bin_index<explicit_ff>(this->value.ff, v7.value.ff),
            xyzff::ff_bin_index<explicit_ff>(this->value.ff, v8.value.ff)
        }
    );
}

template<bool vbw, bool explicit_ff>
inline ausaxs::hist::detail::xyzff::OctoEvaluatedResultRounded ausaxs::hist::detail::CompactCoordinatesXYZFF<vbw, explicit_ff>::evaluate_rounded_scalar(
    const CompactCoordinatesXYZFF& v1, const CompactCoordinatesXYZFF& v2, const CompactCoordinatesXYZFF& v3, const CompactCoordinatesXYZFF& v4, 
    const CompactCoordinatesXYZFF& v5, const CompactCoordinatesXYZFF& v6, const CompactCoordinatesXYZFF& v7, const CompactCoordinatesXYZFF& v8
) const noexcept {
    int32_t dx1 = std::round(get_inv_width()*std::sqrt(squared_dot_product(this->data.data(), v1.data.data())));
    int32_t dx2 = std::round(get_inv_width()*std::sqrt(squared_dot_product(this->data.data(), v2.data.data())));
    int32_t dx3 = std::round(get_inv_width()*std::sqrt(squared_dot_product(this->data.data(), v3.data.data())));
    int32_t dx4 = std::round(get_inv_width()*std::sqrt(squared_dot_product(this->data.data(), v4.data.data())));
    int32_t dx5 = std::round(get_inv_width()*std::sqrt(squared_dot_product(this->data.data(), v5.data.data())));
    int32_t dx6 = std::round(get_inv_width()*std::sqrt(squared_dot_product(this->data.data(), v6.data.data())));
    int32_t dx7 = std::round(get_inv_width()*std::sqrt(squared_dot_product(this->data.data(), v7.data.data())));
    int32_t dx8 = std::round(get_inv_width()*std::sqrt(squared_dot_product(this->data.data(), v8.data.data())));
    return xyzff::OctoEvaluatedResultRounded(
        std::array<int32_t, 8>{dx1, dx2, dx3, dx4, dx5, dx6, dx7, dx8},
        std::array<int32_t, 8>{
            xyzff::ff_bin_index<explicit_ff>(this->value.ff, v1.value.ff),
            xyzff::ff_bin_index<explicit_ff>(this->value.ff, v2.value.ff),
            xyzff::ff_bin_index<explicit_ff>(this->value.ff, v3.value.ff),
            xyzff::ff_bin_index<explicit_ff>(this->value.ff, v4.value.ff),
            xyzff::ff_bin_index<explicit_ff>(this->value.ff, v5.value.ff),
            xyzff::ff_bin_index<explicit_ff>(this->value.ff, v6.value.ff),
            xyzff::ff_bin_index<explicit_ff>(this->value.ff, v7.value.ff),
            xyzff::ff_bin_index<explicit_ff>(this->value.ff, v8.value.ff)
        }
    );
}

#if defined AUSAXS_USE_SSE2
    #include <nmmintrin.h>
    namespace ausaxs::hist::detail::xyzff {
        template<bool vbw, bool explicit_ff>
        inline static __m128 ff_bin_index(
            const CompactCoordinatesXYZFF<vbw, explicit_ff>& v, 
            const CompactCoordinatesXYZFF<vbw, explicit_ff>& v1, const CompactCoordinatesXYZFF<vbw, explicit_ff>& v2, const CompactCoordinatesXYZFF<vbw, explicit_ff>& v3, const CompactCoordinatesXYZFF<vbw, explicit_ff>& v4) noexcept {
            __m128 mul_fac;
            if constexpr (explicit_ff) {
                mul_fac = _mm_set_ps1(form_factor::get_count_without_excluded_volume());
            } else {
                mul_fac = _mm_set_ps1(form_factor::get_count());
            }
            __m128 ff1 = _mm_set_ps1(v.value.ff);
            __m128 ff2 = _mm_set_ps(v4.value.ff, v3.value.ff, v2.value.ff, v1.value.ff);
            __m128 ff1_scaled = _mm_mul_ps(ff1, mul_fac);
            return _mm_add_ps(ff2, ff1_scaled);
        }
    }

    template<bool vbw, bool explicit_ff>
    inline ausaxs::hist::detail::xyzff::EvaluatedResult ausaxs::hist::detail::CompactCoordinatesXYZFF<vbw, explicit_ff>::evaluate_sse(
        const CompactCoordinatesXYZFF& other
    ) const noexcept {
        __m128 dist2 = squared_dot_product(this->data.data(), other.data.data(), OutputControl::ALL);
        __m128 dist_sqrt = _mm_sqrt_ps(dist2);
        float dist = _mm_cvtss_f32(dist_sqrt);
        int32_t dist_bin = std::round(get_inv_width() * dist);
        int32_t ff_bin = xyzff::ff_bin_index<explicit_ff>(this->value.ff, other.value.ff);
        return xyzff::EvaluatedResult(dist, dist_bin, ff_bin);
    }

    template<bool vbw, bool explicit_ff>
    inline ausaxs::hist::detail::xyzff::EvaluatedResultRounded ausaxs::hist::detail::CompactCoordinatesXYZFF<vbw, explicit_ff>::evaluate_rounded_sse(
        const CompactCoordinatesXYZFF& other
    ) const noexcept {
        __m128 dist2 = squared_dot_product(this->data.data(), other.data.data(), OutputControl::ALL);
        __m128 dist_sqrt = _mm_sqrt_ps(dist2);
        int32_t dist_bin = std::round(get_inv_width()*_mm_cvtss_f32(dist_sqrt));
        int32_t ff_bin = xyzff::ff_bin_index<explicit_ff>(this->value.ff, other.value.ff);
        return xyzff::EvaluatedResultRounded(dist_bin, ff_bin);
    }

    template<bool vbw, bool explicit_ff>
    inline ausaxs::hist::detail::xyzff::QuadEvaluatedResult ausaxs::hist::detail::CompactCoordinatesXYZFF<vbw, explicit_ff>::evaluate_sse(
        const CompactCoordinatesXYZFF& v1, const CompactCoordinatesXYZFF& v2, const CompactCoordinatesXYZFF& v3, const CompactCoordinatesXYZFF& v4
    ) const noexcept {
        __m128 sv = _mm_load_ps(this->data.data());
        __m128 d1 = _mm_sub_ps(sv, _mm_load_ps(v1.data.data()));
        __m128 d2 = _mm_sub_ps(sv, _mm_load_ps(v2.data.data()));
        __m128 d3 = _mm_sub_ps(sv, _mm_load_ps(v3.data.data()));
        __m128 d4 = _mm_sub_ps(sv, _mm_load_ps(v4.data.data()));
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
        __m128 ff_bins_f = xyzff::ff_bin_index<vbw, explicit_ff>(*this, v1, v2, v3, v4);
        __m128i ff_bins = _mm_cvtps_epi32(ff_bins_f);

        xyzff::QuadEvaluatedResult result;
        _mm_store_ps(reinterpret_cast<float*>(result.distances.data()), dist_sqrt);
        _mm_store_si128(reinterpret_cast<__m128i*>(result.distance_bins.data()), dist_bin);
        _mm_store_si128(reinterpret_cast<__m128i*>(result.ff_bins.data()), ff_bins);
        return result;
    }

    template<bool vbw, bool explicit_ff>
    inline ausaxs::hist::detail::xyzff::QuadEvaluatedResultRounded ausaxs::hist::detail::CompactCoordinatesXYZFF<vbw, explicit_ff>::evaluate_rounded_sse(
        const CompactCoordinatesXYZFF& v1, const CompactCoordinatesXYZFF& v2, const CompactCoordinatesXYZFF& v3, const CompactCoordinatesXYZFF& v4
    ) const noexcept {
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
        __m128 ff_bins_f = xyzff::ff_bin_index<vbw, explicit_ff>(*this, v1, v2, v3, v4);
        __m128i ff_bins = _mm_cvtps_epi32(ff_bins_f);

        xyzff::QuadEvaluatedResultRounded result;
        _mm_store_si128(reinterpret_cast<__m128i*>(result.distances.data()), dist_bin);
        _mm_store_si128(reinterpret_cast<__m128i*>(result.ff_bins.data()), ff_bins);
        return result;
    }

    template<bool vbw, bool explicit_ff>
    inline ausaxs::hist::detail::xyzff::OctoEvaluatedResult ausaxs::hist::detail::CompactCoordinatesXYZFF<vbw, explicit_ff>::evaluate_sse(
        const CompactCoordinatesXYZFF& v1, const CompactCoordinatesXYZFF& v2, const CompactCoordinatesXYZFF& v3, const CompactCoordinatesXYZFF& v4, 
        const CompactCoordinatesXYZFF& v5, const CompactCoordinatesXYZFF& v6, const CompactCoordinatesXYZFF& v7, const CompactCoordinatesXYZFF& v8
    ) const noexcept {
        __m128 sv = _mm_load_ps(this->data.data());
        xyzff::OctoEvaluatedResult result;
        {   // first four
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
            __m128 ff_bins_f = xyzff::ff_bin_index<vbw, explicit_ff>(*this, v1, v2, v3, v4);
            __m128i ff_bins = _mm_cvtps_epi32(ff_bins_f);

            _mm_store_ps(reinterpret_cast<float*>(result.distances.data()), dist_sqrt);
            _mm_store_si128(reinterpret_cast<__m128i*>(result.distance_bins.data()), dist_bin);
            _mm_store_si128(reinterpret_cast<__m128i*>(result.ff_bins.data()), ff_bins);
        }
        {   // last four
            __m128 d1 = _mm_sub_ps(sv, _mm_load_ps(v5.data.data()));
            __m128 d2 = _mm_sub_ps(sv, _mm_load_ps(v6.data.data()));
            __m128 d3 = _mm_sub_ps(sv, _mm_load_ps(v7.data.data()));
            __m128 d4 = _mm_sub_ps(sv, _mm_load_ps(v8.data.data()));
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
            __m128 ff_bins_f = xyzff::ff_bin_index<vbw, explicit_ff>(*this, v5, v6, v7, v8);
            __m128i ff_bins = _mm_cvtps_epi32(ff_bins_f);

            _mm_store_ps(reinterpret_cast<float*>(result.distances.data()+4), dist_sqrt);
            _mm_store_si128(reinterpret_cast<__m128i*>(result.distance_bins.data()+4), dist_bin);
            _mm_store_si128(reinterpret_cast<__m128i*>(result.ff_bins.data()+4), ff_bins);
        }
        return result;
    }

    template<bool vbw, bool explicit_ff>
    inline ausaxs::hist::detail::xyzff::OctoEvaluatedResultRounded ausaxs::hist::detail::CompactCoordinatesXYZFF<vbw, explicit_ff>::evaluate_rounded_sse(
        const CompactCoordinatesXYZFF& v1, const CompactCoordinatesXYZFF& v2, const CompactCoordinatesXYZFF& v3, const CompactCoordinatesXYZFF& v4, 
        const CompactCoordinatesXYZFF& v5, const CompactCoordinatesXYZFF& v6, const CompactCoordinatesXYZFF& v7, const CompactCoordinatesXYZFF& v8
    ) const noexcept {
        __m128 sv = _mm_load_ps(this->data.data());
        xyzff::OctoEvaluatedResultRounded result;
        {   // first four
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
            __m128 ff_bins_f = xyzff::ff_bin_index<vbw, explicit_ff>(*this, v1, v2, v3, v4);
            __m128i ff_bins = _mm_cvtps_epi32(ff_bins_f);

            _mm_store_si128(reinterpret_cast<__m128i*>(result.distances.data()), dist_bin);
            _mm_store_si128(reinterpret_cast<__m128i*>(result.ff_bins.data()), ff_bins);
        }
        {   // last four
            __m128 d1 = _mm_sub_ps(sv, _mm_load_ps(v5.data.data()));
            __m128 d2 = _mm_sub_ps(sv, _mm_load_ps(v6.data.data()));
            __m128 d3 = _mm_sub_ps(sv, _mm_load_ps(v7.data.data()));
            __m128 d4 = _mm_sub_ps(sv, _mm_load_ps(v8.data.data()));
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
            __m128 ff_bins_f = xyzff::ff_bin_index<vbw, explicit_ff>(*this, v5, v6, v7, v8);
            __m128i ff_bins = _mm_cvtps_epi32(ff_bins_f);

            _mm_store_si128(reinterpret_cast<__m128i*>(result.distances.data()+4), dist_bin);
            _mm_store_si128(reinterpret_cast<__m128i*>(result.ff_bins.data()+4), ff_bins);
        }
        return result;
    }
#endif

#if defined AUSAXS_USE_AVX2
    #include <immintrin.h>
    namespace ausaxs::hist::detail::xyzff {
        template<bool vbw, bool explicit_ff>
        inline static __m256 ff_bin_index(
            const CompactCoordinatesXYZFF<vbw, explicit_ff>& v, 
            const CompactCoordinatesXYZFF<vbw, explicit_ff>& v1, const CompactCoordinatesXYZFF<vbw, explicit_ff>& v2, const CompactCoordinatesXYZFF<vbw, explicit_ff>& v3, const CompactCoordinatesXYZFF<vbw, explicit_ff>& v4,
            const CompactCoordinatesXYZFF<vbw, explicit_ff>& v5, const CompactCoordinatesXYZFF<vbw, explicit_ff>& v6, const CompactCoordinatesXYZFF<vbw, explicit_ff>& v7, const CompactCoordinatesXYZFF<vbw, explicit_ff>& v8
        ) noexcept {
            __m256 mul_fac;
            if constexpr (explicit_ff) {
                mul_fac = _mm256_set1_ps(form_factor::get_count_without_excluded_volume());
            } else {
                mul_fac = _mm256_set1_ps(form_factor::get_count());
            }
            __m256 ff1 = _mm256_set1_ps(v.value.ff);
            __m256 ff2 = _mm256_set_ps(v8.value.ff, v7.value.ff, v6.value.ff, v5.value.ff, v4.value.ff, v3.value.ff, v2.value.ff, v1.value.ff);
            __m256 ff1_scaled = _mm256_mul_ps(ff1, mul_fac);
            return _mm256_add_ps(ff2, ff1_scaled);
        }
    }

    template<bool vbw, bool explicit_ff>
    inline ausaxs::hist::detail::xyzff::EvaluatedResult ausaxs::hist::detail::CompactCoordinatesXYZFF<vbw, explicit_ff>::evaluate_avx(
        const CompactCoordinatesXYZFF& other
    ) const noexcept {
        return evaluate_sse(other); // no way to optimize a single evaluation with AVX
    }

    template<bool vbw, bool explicit_ff>
    inline ausaxs::hist::detail::xyzff::EvaluatedResultRounded ausaxs::hist::detail::CompactCoordinatesXYZFF<vbw, explicit_ff>::evaluate_rounded_avx(
        const CompactCoordinatesXYZFF& other
    ) const noexcept {
        return evaluate_rounded_sse(other); // no way to optimize a single evaluation with AVX
    }

    template<bool vbw, bool explicit_ff>
    inline ausaxs::hist::detail::xyzff::QuadEvaluatedResult ausaxs::hist::detail::CompactCoordinatesXYZFF<vbw, explicit_ff>::evaluate_avx(
        const CompactCoordinatesXYZFF& v1, const CompactCoordinatesXYZFF& v2, const CompactCoordinatesXYZFF& v3, const CompactCoordinatesXYZFF& v4
    ) const noexcept {
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
        __m128 ff_bins_f = xyzff::ff_bin_index<vbw, explicit_ff>(*this, v1, v2, v3, v4);
        __m128i ff_bins = _mm_cvtps_epi32(ff_bins_f);

        xyzff::QuadEvaluatedResult result;
        _mm_store_ps(reinterpret_cast<float*>(result.distances.data()), dist_sqrt);
        _mm_store_si128(reinterpret_cast<__m128i*>(result.distance_bins.data()), dist_bin);
        _mm_store_si128(reinterpret_cast<__m128i*>(result.ff_bins.data()), ff_bins);
        return result;
    }

    template<bool vbw, bool explicit_ff>
    inline ausaxs::hist::detail::xyzff::QuadEvaluatedResultRounded ausaxs::hist::detail::CompactCoordinatesXYZFF<vbw, explicit_ff>::evaluate_rounded_avx(
        const CompactCoordinatesXYZFF& v1, const CompactCoordinatesXYZFF& v2, const CompactCoordinatesXYZFF& v3, const CompactCoordinatesXYZFF& v4
    ) const noexcept {
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
        __m128 ff_bins_f = xyzff::ff_bin_index<vbw, explicit_ff>(*this, v1, v2, v3, v4);
        __m128i ff_bins = _mm_cvtps_epi32(ff_bins_f);

        xyzff::QuadEvaluatedResultRounded result;
        _mm_store_si128(reinterpret_cast<__m128i*>(result.distances.data()), dist_bin);
        _mm_store_si128(reinterpret_cast<__m128i*>(result.ff_bins.data()), ff_bins);
        return result;
    }

    template<bool vbw, bool explicit_ff>
    inline ausaxs::hist::detail::xyzff::OctoEvaluatedResult ausaxs::hist::detail::CompactCoordinatesXYZFF<vbw, explicit_ff>::evaluate_avx(
        const CompactCoordinatesXYZFF& v1, const CompactCoordinatesXYZFF& v2, const CompactCoordinatesXYZFF& v3, const CompactCoordinatesXYZFF& v4, 
        const CompactCoordinatesXYZFF& v5, const CompactCoordinatesXYZFF& v6, const CompactCoordinatesXYZFF& v7, const CompactCoordinatesXYZFF& v8
    ) const noexcept {
        // 4x 256-bit loads from contiguous memory (v1..v8 must be adjacent)
        const float* base = v1.data.data();
        __m256 v12 = _mm256_loadu_ps(base);
        __m256 v34 = _mm256_loadu_ps(base + 8);
        __m256 v56 = _mm256_loadu_ps(base + 16);
        __m256 v78 = _mm256_loadu_ps(base + 24);

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

        xyzff::OctoEvaluatedResult result;
        _mm256_store_ps(reinterpret_cast<float*>(result.distances.data()), dist_sqrt);
        _mm256_store_si256(reinterpret_cast<__m256i*>(result.distance_bins.data()), dist_bin);
        _mm256_store_si256(reinterpret_cast<__m256i*>(result.ff_bins.data()), ff_bins);
        return result;
    }

    template<bool vbw, bool explicit_ff>
    inline ausaxs::hist::detail::xyzff::OctoEvaluatedResultRounded ausaxs::hist::detail::CompactCoordinatesXYZFF<vbw, explicit_ff>::evaluate_rounded_avx(
        const CompactCoordinatesXYZFF& v1, const CompactCoordinatesXYZFF& v2, const CompactCoordinatesXYZFF& v3, const CompactCoordinatesXYZFF& v4, 
        const CompactCoordinatesXYZFF& v5, const CompactCoordinatesXYZFF& v6, const CompactCoordinatesXYZFF& v7, const CompactCoordinatesXYZFF& v8
    ) const noexcept {
        const float* base = v1.data.data();
        __m256 v12 = _mm256_loadu_ps(base);
        __m256 v34 = _mm256_loadu_ps(base + 8);
        __m256 v56 = _mm256_loadu_ps(base + 16);
        __m256 v78 = _mm256_loadu_ps(base + 24);

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

        xyzff::OctoEvaluatedResultRounded result;
        _mm256_store_si256(reinterpret_cast<__m256i*>(result.distances.data()), dist_bin);
        _mm256_store_si256(reinterpret_cast<__m256i*>(result.ff_bins.data()), ff_bins);
        return result;
    }
#endif

#if defined AUSAXS_USE_AVX512
    #include <immintrin.h>

    template<bool vbw, bool explicit_ff>
    inline ausaxs::hist::detail::xyzff::OctoEvaluatedResult ausaxs::hist::detail::CompactCoordinatesXYZFF<vbw, explicit_ff>::evaluate_avx512(
        const CompactCoordinatesXYZFF& v1, const CompactCoordinatesXYZFF& v2, const CompactCoordinatesXYZFF& v3, const CompactCoordinatesXYZFF& v4, 
        const CompactCoordinatesXYZFF& v5, const CompactCoordinatesXYZFF& v6, const CompactCoordinatesXYZFF& v7, const CompactCoordinatesXYZFF& v8
    ) const noexcept {
        const float* base = v1.data.data();
        __m512 v1234 = _mm512_loadu_ps(base);
        __m512 v5678 = _mm512_loadu_ps(base + 16);

        // extract ff values (int32 at position 3 in each atom)
        const __m512i gather_w = _mm512_setr_epi32(3, 7, 11, 15, 19, 23, 27, 31, 0, 0, 0, 0, 0, 0, 0, 0);
        __m256 ff_raw = _mm512_castps512_ps256(_mm512_permutex2var_ps(v1234, gather_w, v5678));
        __m256 ff_float = _mm256_cvtepi32_ps(_mm256_castps_si256(ff_raw));
        __m256 mul_fac;
        if constexpr (explicit_ff) {
            mul_fac = _mm256_set1_ps(form_factor::get_count_without_excluded_volume());
        } else {
            mul_fac = _mm256_set1_ps(form_factor::get_count());
        }
        __m256i ff_bins = _mm256_cvtps_epi32(_mm256_add_ps(ff_float, _mm256_mul_ps(_mm256_set1_ps(this->value.ff), mul_fac)));

        __m512 svv = _mm512_broadcast_f32x4(_mm_load_ps(this->data.data()));
        __m512 d1234 = _mm512_sub_ps(svv, v1234);
        __m512 d5678 = _mm512_sub_ps(svv, v5678);

        // cross-lane gather of dx, dy, dz from all 8 atoms, then FMA dist²
        const __m512i gather_x = _mm512_setr_epi32(0, 4, 8, 12, 16, 20, 24, 28, 0, 0, 0, 0, 0, 0, 0, 0);
        const __m512i gather_y = _mm512_setr_epi32(1, 5, 9, 13, 17, 21, 25, 29, 0, 0, 0, 0, 0, 0, 0, 0);
        const __m512i gather_z = _mm512_setr_epi32(2, 6, 10, 14, 18, 22, 26, 30, 0, 0, 0, 0, 0, 0, 0, 0);
        __m256 dx = _mm512_castps512_ps256(_mm512_permutex2var_ps(d1234, gather_x, d5678));
        __m256 dy = _mm512_castps512_ps256(_mm512_permutex2var_ps(d1234, gather_y, d5678));
        __m256 dz = _mm512_castps512_ps256(_mm512_permutex2var_ps(d1234, gather_z, d5678));

        __m256 dist2 = _mm256_mul_ps(dx, dx);
        dist2 = _mm256_fmadd_ps(dy, dy, dist2);
        dist2 = _mm256_fmadd_ps(dz, dz, dist2);

        __m256 dist_sqrt = _mm256_sqrt_ps(dist2);
        __m256i dist_bin = _mm256_cvtps_epi32(_mm256_mul_ps(dist_sqrt, _mm256_set1_ps(get_inv_width())));

        xyzff::OctoEvaluatedResult result;
        _mm256_store_ps(reinterpret_cast<float*>(result.distances.data()), dist_sqrt);
        _mm256_store_si256(reinterpret_cast<__m256i*>(result.distance_bins.data()), dist_bin);
        _mm256_store_si256(reinterpret_cast<__m256i*>(result.ff_bins.data()), ff_bins);
        return result;
    }

    template<bool vbw, bool explicit_ff>
    inline ausaxs::hist::detail::xyzff::OctoEvaluatedResultRounded ausaxs::hist::detail::CompactCoordinatesXYZFF<vbw, explicit_ff>::evaluate_rounded_avx512(
        const CompactCoordinatesXYZFF& v1, const CompactCoordinatesXYZFF& v2, const CompactCoordinatesXYZFF& v3, const CompactCoordinatesXYZFF& v4, 
        const CompactCoordinatesXYZFF& v5, const CompactCoordinatesXYZFF& v6, const CompactCoordinatesXYZFF& v7, const CompactCoordinatesXYZFF& v8
    ) const noexcept {
        const float* base = v1.data.data();
        __m512 v1234 = _mm512_loadu_ps(base);
        __m512 v5678 = _mm512_loadu_ps(base + 16);

        const __m512i gather_w = _mm512_setr_epi32(3, 7, 11, 15, 19, 23, 27, 31, 0, 0, 0, 0, 0, 0, 0, 0);
        __m256 ff_raw = _mm512_castps512_ps256(_mm512_permutex2var_ps(v1234, gather_w, v5678));
        __m256 ff_float = _mm256_cvtepi32_ps(_mm256_castps_si256(ff_raw));
        __m256 mul_fac;
        if constexpr (explicit_ff) {
            mul_fac = _mm256_set1_ps(form_factor::get_count_without_excluded_volume());
        } else {
            mul_fac = _mm256_set1_ps(form_factor::get_count());
        }
        __m256i ff_bins = _mm256_cvtps_epi32(_mm256_add_ps(ff_float, _mm256_mul_ps(_mm256_set1_ps(this->value.ff), mul_fac)));

        __m512 svv = _mm512_broadcast_f32x4(_mm_load_ps(this->data.data()));
        __m512 d1234 = _mm512_sub_ps(svv, v1234);
        __m512 d5678 = _mm512_sub_ps(svv, v5678);

        const __m512i gather_x = _mm512_setr_epi32(0, 4, 8, 12, 16, 20, 24, 28, 0, 0, 0, 0, 0, 0, 0, 0);
        const __m512i gather_y = _mm512_setr_epi32(1, 5, 9, 13, 17, 21, 25, 29, 0, 0, 0, 0, 0, 0, 0, 0);
        const __m512i gather_z = _mm512_setr_epi32(2, 6, 10, 14, 18, 22, 26, 30, 0, 0, 0, 0, 0, 0, 0, 0);
        __m256 dx = _mm512_castps512_ps256(_mm512_permutex2var_ps(d1234, gather_x, d5678));
        __m256 dy = _mm512_castps512_ps256(_mm512_permutex2var_ps(d1234, gather_y, d5678));
        __m256 dz = _mm512_castps512_ps256(_mm512_permutex2var_ps(d1234, gather_z, d5678));

        __m256 dist2 = _mm256_mul_ps(dx, dx);
        dist2 = _mm256_fmadd_ps(dy, dy, dist2);
        dist2 = _mm256_fmadd_ps(dz, dz, dist2);

        __m256i dist_bin = _mm256_cvtps_epi32(_mm256_mul_ps(
            _mm256_sqrt_ps(dist2),
            _mm256_set1_ps(get_inv_width())));

        xyzff::OctoEvaluatedResultRounded result;
        _mm256_store_si256(reinterpret_cast<__m256i*>(result.distances.data()), dist_bin);
        _mm256_store_si256(reinterpret_cast<__m256i*>(result.ff_bins.data()), ff_bins);
        return result;
    }
#endif