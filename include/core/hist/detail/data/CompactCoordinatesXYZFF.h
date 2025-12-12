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

// AVX implies SSE4.1 and SSE2, but the MSVC compiler doesn't seem to define the latter two
#if defined __AVX__
    #if !defined __SSE2__
        #define __SSE2__
    #endif
    #if !defined __SSE4_1__
        #define __SSE4_1__
    #endif
#endif

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
    struct OctoEvaluatedResult {
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
    struct OctoEvaluatedResultRounded {
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

            #if defined __SSE2__
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

            #if defined __AVX__
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

// undefine SSE2 & AVX for MacOS
#if defined __APPLE__
    #undef __AVX__
    #undef __SSE2__
    #undef __SSE4_2__
#endif

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
    #if defined __SSE2__
        return evaluate_sse(other);
    #else
        return evaluate_scalar(other);
    #endif
}

template<bool vbw, bool explicit_ff>
inline ausaxs::hist::detail::xyzff::EvaluatedResultRounded ausaxs::hist::detail::CompactCoordinatesXYZFF<vbw, explicit_ff>::evaluate_rounded(const CompactCoordinatesXYZFF& other) const  noexcept{
    #if defined __SSE2__
        return evaluate_rounded_sse(other);
    #else
        return evaluate_rounded_scalar(other);
    #endif
}

template<bool vbw, bool explicit_ff>
inline ausaxs::hist::detail::xyzff::QuadEvaluatedResult ausaxs::hist::detail::CompactCoordinatesXYZFF<vbw, explicit_ff>::evaluate(
    const CompactCoordinatesXYZFF& v1, const CompactCoordinatesXYZFF& v2, const CompactCoordinatesXYZFF& v3, const CompactCoordinatesXYZFF& v4
) const noexcept {
    #if defined __AVX__
        return evaluate_avx(v1, v2, v3, v4);
    #elif defined __SSE2__
        return evaluate_sse(v1, v2, v3, v4);
    #else
        return evaluate_scalar(v1, v2, v3, v4);
    #endif
}

template<bool vbw, bool explicit_ff>
inline ausaxs::hist::detail::xyzff::QuadEvaluatedResultRounded ausaxs::hist::detail::CompactCoordinatesXYZFF<vbw, explicit_ff>::evaluate_rounded(
    const CompactCoordinatesXYZFF& v1, const CompactCoordinatesXYZFF& v2, const CompactCoordinatesXYZFF& v3, const CompactCoordinatesXYZFF& v4
) const noexcept {
    #if defined __AVX__
        return evaluate_rounded_avx(v1, v2, v3, v4);
    #elif defined __SSE2__
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
    #if defined __AVX__
        return evaluate_avx(v1, v2, v3, v4, v5, v6, v7, v8);
    #elif defined __SSE2__
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
    #if defined __AVX__
        return evaluate_rounded_avx(v1, v2, v3, v4, v5, v6, v7, v8);
    #elif defined __SSE2__
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

#if defined __SSE2__
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
        __m128 dist2_1 = squared_dot_product(this->data.data(), v1.data.data(), OutputControl::FIRST);
        __m128 dist2_2 = squared_dot_product(this->data.data(), v2.data.data(), OutputControl::SECOND);
        __m128 dist2_3 = squared_dot_product(this->data.data(), v3.data.data(), OutputControl::THIRD);
        __m128 dist2_4 = squared_dot_product(this->data.data(), v4.data.data(), OutputControl::FOURTH);

        // since we masked the dot product to four different sections, we can just add the results to get |Δx1^2|Δx2^2|Δx3^2|Δx4^2|
        __m128 dist2 = _mm_add_ps(dist2_1, dist2_2);
        dist2 = _mm_add_ps(dist2, dist2_3);
        dist2 = _mm_add_ps(dist2, dist2_4);

        __m128 dist_sqrt = _mm_sqrt_ps(dist2);
        __m128 dist_binf = _mm_mul_ps(dist_sqrt, _mm_set_ps1(get_inv_width())); // multiply by the inverse width
        __m128i dist_bin = _mm_cvtps_epi32(dist_binf);                          // convert float distances to int bins
        __m128 ff_bins_f = xyzff::ff_bin_index<vbw, explicit_ff>(*this, v1, v2, v3, v4);
        __m128i ff_bins = _mm_cvtps_epi32(ff_bins_f);                           // convert float ff_bins to int

        xyzff::QuadEvaluatedResult result;
        _mm_store_ps(reinterpret_cast<float*>(result.distances.data()), dist_sqrt);             // store distances
        _mm_store_si128(reinterpret_cast<__m128i*>(result.distance_bins.data()), dist_bin);    // store distance bins
        _mm_store_si128(reinterpret_cast<__m128i*>(result.ff_bins.data()), ff_bins);           // store ff_bins
        return result;
    }

    template<bool vbw, bool explicit_ff>
    inline ausaxs::hist::detail::xyzff::QuadEvaluatedResultRounded ausaxs::hist::detail::CompactCoordinatesXYZFF<vbw, explicit_ff>::evaluate_rounded_sse(
        const CompactCoordinatesXYZFF& v1, const CompactCoordinatesXYZFF& v2, const CompactCoordinatesXYZFF& v3, const CompactCoordinatesXYZFF& v4
    ) const noexcept {
        __m128 dist2_1 = squared_dot_product(this->data.data(), v1.data.data(), OutputControl::FIRST);
        __m128 dist2_2 = squared_dot_product(this->data.data(), v2.data.data(), OutputControl::SECOND);
        __m128 dist2_3 = squared_dot_product(this->data.data(), v3.data.data(), OutputControl::THIRD);
        __m128 dist2_4 = squared_dot_product(this->data.data(), v4.data.data(), OutputControl::FOURTH);

        // since we masked the dot product to four different sections, we can just add the results to get |Δx1^2|Δx2^2|Δx3^2|Δx4^2|
        __m128 dist2 = _mm_add_ps(dist2_1, dist2_2);
        dist2 = _mm_add_ps(dist2, dist2_3);
        dist2 = _mm_add_ps(dist2, dist2_4);

        __m128 dist_sqrt = _mm_sqrt_ps(dist2);
        __m128 dist_binf = _mm_mul_ps(dist_sqrt, _mm_set_ps1(get_inv_width())); // multiply by the inverse width
        __m128i dist_bin = _mm_cvtps_epi32(dist_binf);                          // convert float distances to int bins
        __m128 ff_bins_f = xyzff::ff_bin_index<vbw, explicit_ff>(*this, v1, v2, v3, v4);
        __m128i ff_bins = _mm_cvtps_epi32(ff_bins_f);                           // convert float ff_bins to int

        xyzff::QuadEvaluatedResultRounded result;
        _mm_store_si128(reinterpret_cast<__m128i*>(result.distances.data()), dist_bin); // store distance bins
        _mm_store_si128(reinterpret_cast<__m128i*>(result.ff_bins.data()), ff_bins);    // store ff_bins
        return result;
    }

    template<bool vbw, bool explicit_ff>
    inline ausaxs::hist::detail::xyzff::OctoEvaluatedResult ausaxs::hist::detail::CompactCoordinatesXYZFF<vbw, explicit_ff>::evaluate_sse(
        const CompactCoordinatesXYZFF& v1, const CompactCoordinatesXYZFF& v2, const CompactCoordinatesXYZFF& v3, const CompactCoordinatesXYZFF& v4, 
        const CompactCoordinatesXYZFF& v5, const CompactCoordinatesXYZFF& v6, const CompactCoordinatesXYZFF& v7, const CompactCoordinatesXYZFF& v8
    ) const noexcept {
        xyzff::OctoEvaluatedResult result;
        __m128 inv_width = _mm_set_ps1(get_inv_width());
        {   // first four
            __m128 dist2_1 = squared_dot_product(this->data.data(), v1.data.data(), OutputControl::FIRST);
            __m128 dist2_2 = squared_dot_product(this->data.data(), v2.data.data(), OutputControl::SECOND);
            __m128 dist2_3 = squared_dot_product(this->data.data(), v3.data.data(), OutputControl::THIRD);
            __m128 dist2_4 = squared_dot_product(this->data.data(), v4.data.data(), OutputControl::FOURTH);

            __m128 dist2 = _mm_add_ps(dist2_1, dist2_2);
            dist2 = _mm_add_ps(dist2, dist2_3);
            dist2 = _mm_add_ps(dist2, dist2_4);
            __m128 dist_sqrt = _mm_sqrt_ps(dist2);
            __m128 dist_binf = _mm_mul_ps(dist_sqrt, inv_width);              // multiply by the inverse width
            __m128i dist_bin = _mm_cvtps_epi32(dist_binf);                    // convert float distances to int bins
            __m128 ff_bins_f = xyzff::ff_bin_index<vbw, explicit_ff>(*this, v1, v2, v3, v4);
            __m128i ff_bins = _mm_cvtps_epi32(ff_bins_f);

            _mm_store_ps(reinterpret_cast<float*>(result.distances.data()), dist_sqrt);
            _mm_store_si128(reinterpret_cast<__m128i*>(result.distance_bins.data()), dist_bin);
            _mm_store_si128(reinterpret_cast<__m128i*>(result.ff_bins.data()), ff_bins);
        }
        {   // last four
            __m128 dist2_1 = squared_dot_product(this->data.data(), v5.data.data(), OutputControl::FIRST);
            __m128 dist2_2 = squared_dot_product(this->data.data(), v6.data.data(), OutputControl::SECOND);
            __m128 dist2_3 = squared_dot_product(this->data.data(), v7.data.data(), OutputControl::THIRD);
            __m128 dist2_4 = squared_dot_product(this->data.data(), v8.data.data(), OutputControl::FOURTH);

            __m128 dist2 = _mm_add_ps(dist2_1, dist2_2);
            dist2 = _mm_add_ps(dist2, dist2_3);
            dist2 = _mm_add_ps(dist2, dist2_4);
            __m128 dist_sqrt = _mm_sqrt_ps(dist2);
            __m128 dist_binf = _mm_mul_ps(dist_sqrt, inv_width);              // multiply by the inverse width
            __m128i dist_bin = _mm_cvtps_epi32(dist_binf);                    // convert float distances to int bins
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
        xyzff::OctoEvaluatedResultRounded result;
        {   // first four
            __m128 dist2_1 = squared_dot_product(this->data.data(), v1.data.data(), OutputControl::FIRST);
            __m128 dist2_2 = squared_dot_product(this->data.data(), v2.data.data(), OutputControl::SECOND);
            __m128 dist2_3 = squared_dot_product(this->data.data(), v3.data.data(), OutputControl::THIRD);
            __m128 dist2_4 = squared_dot_product(this->data.data(), v4.data.data(), OutputControl::FOURTH);

            __m128 dist2 = _mm_add_ps(dist2_1, dist2_2);
            dist2 = _mm_add_ps(dist2, dist2_3);
            dist2 = _mm_add_ps(dist2, dist2_4);
            __m128 dist_sqrt = _mm_sqrt_ps(dist2);
            __m128 dist_binf = _mm_mul_ps(dist_sqrt, _mm_set_ps1(get_inv_width()));   // multiply by the inverse width
            __m128i dist_bin = _mm_cvtps_epi32(dist_binf);                            // convert float distances to int bins
            __m128 ff_bins_f = xyzff::ff_bin_index<explicit_ff>(*this, v1, v2, v3, v4);
            __m128i ff_bins = _mm_cvtps_epi32(ff_bins_f);

            _mm_store_si128(reinterpret_cast<__m128i*>(result.distances.data()), dist_bin);
            _mm_store_si128(reinterpret_cast<__m128i*>(result.ff_bins.data()), ff_bins);
        }
        {   // last four
            __m128 dist2_1 = squared_dot_product(this->data.data(), v5.data.data(), OutputControl::FIRST);
            __m128 dist2_2 = squared_dot_product(this->data.data(), v6.data.data(), OutputControl::SECOND);
            __m128 dist2_3 = squared_dot_product(this->data.data(), v7.data.data(), OutputControl::THIRD);
            __m128 dist2_4 = squared_dot_product(this->data.data(), v8.data.data(), OutputControl::FOURTH);

            __m128 dist2 = _mm_add_ps(dist2_1, dist2_2);
            dist2 = _mm_add_ps(dist2, dist2_3);
            dist2 = _mm_add_ps(dist2, dist2_4);
            __m128 dist_sqrt = _mm_sqrt_ps(dist2);
            __m128 dist_binf = _mm_mul_ps(dist_sqrt, _mm_set_ps1(get_inv_width()));   // multiply by the inverse width
            __m128i dist_bin = _mm_cvtps_epi32(dist_binf);                            // convert float distances to int bins
            __m128 ff_bins_f = xyzff::ff_bin_index<vbw, explicit_ff>(*this, v5, v6, v7, v8);
            __m128i ff_bins = _mm_cvtps_epi32(ff_bins_f);

            _mm_store_si128(reinterpret_cast<__m128i*>(result.distances.data()+4), dist_bin);
            _mm_store_si128(reinterpret_cast<__m128i*>(result.ff_bins.data()+4), ff_bins);
        }
        return result;
    }
#endif

#if defined __AVX__
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
        __m256 dist2_1 = squared_dot_product(this->data.data(), v1.data.data(), v3.data.data(), OutputControl::FIRST); // |Δx1^2|0    |0    |0    |Δx3^2|0    |0    |0    |
        __m256 dist2_2 = squared_dot_product(this->data.data(), v2.data.data(), v4.data.data(), OutputControl::SECOND);// |0    |Δx2^2|0    |0    |0    |Δx4^2|0    |0    |

        __m256 dist2_256 = _mm256_add_ps(dist2_1, dist2_2);                     // |Δx1^2|Δx2^2|0    |0    |Δx3^2|Δx4^2|0    |0    |
        __m256 dist2_256_shuffled = _mm256_permute_ps(dist2_256, 0b01001110);   // |0    |0    |Δx3^2|Δx4^2|0    |0    |Δx1^2|Δx2^2|
        __m128 dist2_128_lower = _mm256_extractf128_ps(dist2_256, 0);           // |Δx1^2|Δx2^2|0    |0    |
        __m128 dist2_128_upper = _mm256_extractf128_ps(dist2_256_shuffled, 1);  // |0    |0    |Δx3^2|Δx4^2|
        __m128 dist2 = _mm_add_ps(dist2_128_lower, dist2_128_upper);            // |Δx1^2|Δx2^2|Δx3^2|Δx4^2|
        __m128 dist_sqrt = _mm_sqrt_ps(dist2);
        __m128 dist_binf = _mm_mul_ps(dist_sqrt, _mm_set_ps1(get_inv_width())); // multiply by the inverse width
        __m128i dist_bin = _mm_cvtps_epi32(dist_binf);                          // convert float distances to int bins
        __m128 ff_bins_f = xyzff::ff_bin_index<vbw, explicit_ff>(*this, v1, v2, v3, v4);
        __m128i ff_bins = _mm_cvtps_epi32(ff_bins_f);                           // convert float ff_bins to int

        xyzff::QuadEvaluatedResult result;
        _mm_store_ps(reinterpret_cast<float*>(result.distances.data()), dist_sqrt);             // store distances
        _mm_store_si128(reinterpret_cast<__m128i*>(result.distance_bins.data()), dist_bin);    // store distance bins
        _mm_store_si128(reinterpret_cast<__m128i*>(result.ff_bins.data()), ff_bins);           // store ff_bins
        return result;
    }

    template<bool vbw, bool explicit_ff>
    inline ausaxs::hist::detail::xyzff::QuadEvaluatedResultRounded ausaxs::hist::detail::CompactCoordinatesXYZFF<vbw, explicit_ff>::evaluate_rounded_avx(
        const CompactCoordinatesXYZFF& v1, const CompactCoordinatesXYZFF& v2, const CompactCoordinatesXYZFF& v3, const CompactCoordinatesXYZFF& v4
    ) const noexcept {
        __m256 dist2_1 = squared_dot_product(this->data.data(), v1.data.data(), v3.data.data(), OutputControl::FIRST); // |Δx1^2|0    |0    |0    |Δx3^2|0    |0    |0    |
        __m256 dist2_2 = squared_dot_product(this->data.data(), v2.data.data(), v4.data.data(), OutputControl::SECOND);// |0    |Δx2^2|0    |0    |0    |Δx4^2|0    |0    |

        __m256 dist2_256 = _mm256_add_ps(dist2_1, dist2_2);                     // |Δx1^2|Δx2^2|0    |0    |Δx3^2|Δx4^2|0    |0    |
        __m256 dist2_256_shuffled = _mm256_permute_ps(dist2_256, 0b01001110);   // |0    |0    |Δx3^2|Δx4^2|0    |0    |Δx1^2|Δx2^2|
        __m128 dist2_128_lower = _mm256_extractf128_ps(dist2_256, 0);           // |Δx1^2|Δx2^2|0    |0    |
        __m128 dist2_128_upper = _mm256_extractf128_ps(dist2_256_shuffled, 1);  // |0    |0    |Δx3^2|Δx4^2|
        __m128 dist2 = _mm_add_ps(dist2_128_lower, dist2_128_upper);            // |Δx1^2|Δx2^2|Δx3^2|Δx4^2|
        __m128 dist_sqrt = _mm_sqrt_ps(dist2);
        __m128 dist_binf = _mm_mul_ps(dist_sqrt, _mm_set_ps1(get_inv_width())); // multiply by the inverse width
        __m128i dist_bin = _mm_cvtps_epi32(dist_binf);                          // convert float distances to int bins
        __m128 ff_bins_f = xyzff::ff_bin_index<vbw, explicit_ff>(*this, v1, v2, v3, v4);
        __m128i ff_bins = _mm_cvtps_epi32(ff_bins_f);                           // convert float ff_bins to int

        xyzff::QuadEvaluatedResultRounded result;
        _mm_store_si128(reinterpret_cast<__m128i*>(result.distances.data()), dist_bin);  // store distance bins
        _mm_store_si128(reinterpret_cast<__m128i*>(result.ff_bins.data()), ff_bins);     // store ff_bins
        return result;
    }

    template<bool vbw, bool explicit_ff>
    inline ausaxs::hist::detail::xyzff::OctoEvaluatedResult ausaxs::hist::detail::CompactCoordinatesXYZFF<vbw, explicit_ff>::evaluate_avx(
        const CompactCoordinatesXYZFF& v1, const CompactCoordinatesXYZFF& v2, const CompactCoordinatesXYZFF& v3, const CompactCoordinatesXYZFF& v4, 
        const CompactCoordinatesXYZFF& v5, const CompactCoordinatesXYZFF& v6, const CompactCoordinatesXYZFF& v7, const CompactCoordinatesXYZFF& v8
    ) const noexcept {
        __m256 dist2_1 = squared_dot_product(this->data.data(), v1.data.data(), v5.data.data(), OutputControl::FIRST); // |Δx1^2|0    |0    |0    |Δx5^2|0    |0    |0    |
        __m256 dist2_2 = squared_dot_product(this->data.data(), v2.data.data(), v6.data.data(), OutputControl::SECOND);// |0    |Δx2^2|0    |0    |0    |Δx6^2|0    |0    |
        __m256 dist2_3 = squared_dot_product(this->data.data(), v3.data.data(), v7.data.data(), OutputControl::THIRD); // |0    |0    |Δx3^2|0    |0    |0    |Δx7^2|0    |
        __m256 dist2_4 = squared_dot_product(this->data.data(), v4.data.data(), v8.data.data(), OutputControl::FOURTH);// |0    |0    |0    |Δx4^2|0    |0    |0    |Δx8^2|

        __m256 dist2_256 = _mm256_add_ps(dist2_1, dist2_2);                                                            // |Δx1^2|Δx2^2|0    |0    |Δx5^2|Δx6^2|0    |0    |
        dist2_256 = _mm256_add_ps(dist2_256, dist2_3);                                                                 // |Δx1^2|Δx2^2|Δx3^2|0    |Δx5^2|Δx6^2|Δx7^2|0    |
        dist2_256 = _mm256_add_ps(dist2_256, dist2_4);                                                                 // |Δx1^2|Δx2^2|Δx3^2|Δx4^2|Δx5^2|Δx6^2|Δx7^2|Δx8^2|
        __m256 dist_sqrt = _mm256_sqrt_ps(dist2_256);
        __m256 dist_binf = _mm256_mul_ps(dist_sqrt, _mm256_set1_ps(get_inv_width())); // multiply by the inverse width
        __m256i dist_bin = _mm256_cvtps_epi32(dist_binf);                             // convert float distances to int bins
        __m256 ff_bins_f = xyzff::ff_bin_index<vbw, explicit_ff>(*this, v1, v2, v3, v4, v5, v6, v7, v8);
        __m256i ff_bins = _mm256_cvtps_epi32(ff_bins_f);                              // convert float ff_bins to int

        xyzff::OctoEvaluatedResult result;
        _mm256_store_ps(reinterpret_cast<float*>(result.distances.data()), dist_sqrt);             // store distances
        _mm256_store_si256(reinterpret_cast<__m256i*>(result.distance_bins.data()), dist_bin);    // store distance bins
        _mm256_store_si256(reinterpret_cast<__m256i*>(result.ff_bins.data()), ff_bins);           // store ff_bins
        return result;
    }

    template<bool vbw, bool explicit_ff>
    inline ausaxs::hist::detail::xyzff::OctoEvaluatedResultRounded ausaxs::hist::detail::CompactCoordinatesXYZFF<vbw, explicit_ff>::evaluate_rounded_avx(
        const CompactCoordinatesXYZFF& v1, const CompactCoordinatesXYZFF& v2, const CompactCoordinatesXYZFF& v3, const CompactCoordinatesXYZFF& v4, 
        const CompactCoordinatesXYZFF& v5, const CompactCoordinatesXYZFF& v6, const CompactCoordinatesXYZFF& v7, const CompactCoordinatesXYZFF& v8
    ) const noexcept {
        __m256 dist2_1 = squared_dot_product(this->data.data(), v1.data.data(), v5.data.data(), OutputControl::FIRST); // |Δx1^2|0    |0    |0    |Δx5^2|0    |0    |0    |
        __m256 dist2_2 = squared_dot_product(this->data.data(), v2.data.data(), v6.data.data(), OutputControl::SECOND);// |0    |Δx2^2|0    |0    |0    |Δx6^2|0    |0    |
        __m256 dist2_3 = squared_dot_product(this->data.data(), v3.data.data(), v7.data.data(), OutputControl::THIRD); // |0    |0    |Δx3^2|0    |0    |0    |Δx7^2|0    |
        __m256 dist2_4 = squared_dot_product(this->data.data(), v4.data.data(), v8.data.data(), OutputControl::FOURTH);// |0    |0    |0    |Δx4^2|0    |0    |0    |Δx8^2|

        __m256 dist2_256 = _mm256_add_ps(dist2_1, dist2_2);                                                            // |Δx1^2|Δx2^2|0    |0    |Δx5^2|Δx6^2|0    |0    |
        dist2_256 = _mm256_add_ps(dist2_256, dist2_3);                                                                 // |Δx1^2|Δx2^2|Δx3^2|0    |Δx5^2|Δx6^2|Δx7^2|0    |
        dist2_256 = _mm256_add_ps(dist2_256, dist2_4);                                                                 // |Δx1^2|Δx2^2|Δx3^2|Δx4^2|Δx5^2|Δx6^2|Δx7^2|Δx8^2|
        __m256 dist_sqrt = _mm256_sqrt_ps(dist2_256);
        __m256 dist_binf = _mm256_mul_ps(dist_sqrt, _mm256_set1_ps(get_inv_width())); // multiply by the inverse width
        __m256i dist_bin = _mm256_cvtps_epi32(dist_binf);                             // convert float distances to int bins
        __m256 ff_bins_f = xyzff::ff_bin_index<vbw, explicit_ff>(*this, v1, v2, v3, v4, v5, v6, v7, v8);
        __m256i ff_bins = _mm256_cvtps_epi32(ff_bins_f);                              // convert float ff_bins to int

        xyzff::OctoEvaluatedResultRounded result;
        _mm256_store_si256(reinterpret_cast<__m256i*>(result.distances.data()), dist_bin);  // store distance bins
        _mm256_store_si256(reinterpret_cast<__m256i*>(result.ff_bins.data()), ff_bins);     // store ff_bins
        return result;
    }
#endif