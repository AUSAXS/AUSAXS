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

namespace ausaxs::hist::detail::xyzw {
    /**
     * @brief Simple structure for storing the results of a distance and weight calculation.
     */
    struct EvaluatedResult {
        EvaluatedResult() noexcept = default;
        EvaluatedResult(float distance, float weight) noexcept : distance(distance), weight(weight) {}
        float distance;         // The raw distance 
        float weight;           // The combined weight
    };

    struct EvaluatedResultRounded {
        EvaluatedResultRounded() noexcept = default;
        EvaluatedResultRounded(int32_t distance_bin, float weight) noexcept : distance(distance_bin), weight(weight) {}
        int32_t distance;       // The distance bin 
        float weight;           // The combined weight
    };

    /**
     * @brief Simple structure for storing the results of four distance and weight calculations.
     *        This is necessary to restructure the memory layout for more efficient 128-bit SIMD instructions.
     */
    struct QuadEvaluatedResult {
        QuadEvaluatedResult() noexcept = default;
        QuadEvaluatedResult(const EvaluatedResult& v1, const EvaluatedResult& v2, const EvaluatedResult& v3, const EvaluatedResult& v4) noexcept 
            : distances{v1.distance, v2.distance, v3.distance, v4.distance}, weights{v1.weight, v2.weight, v3.weight, v4.weight}
        {}
        QuadEvaluatedResult(const std::array<float, 4>& distances, const std::array<float, 4>& weights) noexcept 
            : distances(distances), weights(weights) 
        {}

        std::array<float, 4> distances; // The raw distances
        std::array<float, 4> weights;   // The combined weight
    };

    /**
     * @brief Simple structure for storing the results of four distance and weight calculations.
     *        This is necessary to restructure the memory layout for more efficient 128-bit SIMD instructions.
     */
    struct QuadEvaluatedResultRounded {
        QuadEvaluatedResultRounded() noexcept = default;
        QuadEvaluatedResultRounded(const EvaluatedResultRounded& v1, const EvaluatedResultRounded& v2, const EvaluatedResultRounded& v3, const EvaluatedResultRounded& v4) noexcept 
            : distances{v1.distance, v2.distance, v3.distance, v4.distance}, 
                weights{v1.weight, v2.weight, v3.weight, v4.weight} 
        {}
        QuadEvaluatedResultRounded(const std::array<int32_t, 4>& distances, const std::array<float, 4>& weights) noexcept 
            : distances(distances), weights(weights) 
        {}

        std::array<int32_t, 4> distances;   // The distance bin
        std::array<float, 4> weights;       // The combined weight
    };

    /**
     * @brief Simple structure for storing the results of eight distance and weight calculations.
     *        This is necessary to restructure the memory layout for more efficient 256-bit SIMD instructions.
     */
    struct OctoEvaluatedResult {
        OctoEvaluatedResult() noexcept = default;
        OctoEvaluatedResult(
            const EvaluatedResult& v1, const EvaluatedResult& v2, const EvaluatedResult& v3, const EvaluatedResult& v4, 
            const EvaluatedResult& v5, const EvaluatedResult& v6, const EvaluatedResult& v7, const EvaluatedResult& v8) noexcept 
            : distances{v1.distance, v2.distance, v3.distance, v4.distance, v5.distance, v6.distance, v7.distance, v8.distance},
              weights{v1.weight, v2.weight, v3.weight, v4.weight, v5.weight, v6.weight, v7.weight, v8.weight}
        {}
        OctoEvaluatedResult(const std::array<float, 8>& distances, const std::array<float, 8>& weights) noexcept 
            : distances(distances), weights(weights) 
        {}

        std::array<float, 8> distances; // The distance
        std::array<float, 8> weights;   // The combined weight
    };

    /**
     * @brief Simple structure for storing the results of eight distance and weight calculations.
     *        This is necessary to restructure the memory layout for more efficient 256-bit SIMD instructions.
     */
    struct OctoEvaluatedResultRounded {
        OctoEvaluatedResultRounded() noexcept = default;
        OctoEvaluatedResultRounded(
            const EvaluatedResultRounded& v1, const EvaluatedResultRounded& v2, const EvaluatedResultRounded& v3, const EvaluatedResultRounded& v4, 
            const EvaluatedResultRounded& v5, const EvaluatedResultRounded& v6, const EvaluatedResultRounded& v7, const EvaluatedResultRounded& v8) noexcept 
            : distances{v1.distance, v2.distance, v3.distance, v4.distance, v5.distance, v6.distance, v7.distance, v8.distance}, 
              weights{v1.weight, v2.weight, v3.weight, v4.weight, v5.weight, v6.weight, v7.weight, v8.weight} 
        {}
        OctoEvaluatedResultRounded(const std::array<int32_t, 8>& distances, const std::array<float, 8>& weights) noexcept 
            : distances(distances), weights(weights) 
        {}

        std::array<int32_t, 8> distances;   // The distance bin
        std::array<float, 8> weights;       // The combined weight
    };

    // assert that it is safe to perform memcpy and reinterpret_cast on these structures
    // sizes - EvaluatedResult must be exactly 1 float larger than EvaluatedResultRounded for storing the exact distance
    static_assert(sizeof(EvaluatedResult)            == 8,  "hist::detail::EvaluatedResult is not 12 bytes long");
    static_assert(sizeof(EvaluatedResultRounded)     == 8,  "hist::detail::EvaluatedResultRounded is not 8 bytes long");
    static_assert(sizeof(QuadEvaluatedResult)        == 32, "hist::detail::QuadEvaluatedResult is not 48 bytes long");
    static_assert(sizeof(QuadEvaluatedResultRounded) == 32, "hist::detail::QuadEvaluatedResultRounded is not 32 bytes long");
    static_assert(sizeof(OctoEvaluatedResult)        == 64, "hist::detail::OctoEvaluatedResult is not 96 bytes long");
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

            /**
             * @brief Calculate the @a binned distance and combined weight between this and a single other CompactCoordinatesXYZW.
             */
            xyzw::EvaluatedResultRounded evaluate_rounded(const CompactCoordinatesXYZW& other) const noexcept;

            /**
             * @brief Calculate the distance and combined weight between this and a single other CompactCoordinatesXYZW.
             */
            xyzw::EvaluatedResult evaluate(const CompactCoordinatesXYZW& other) const noexcept;

            /**
             * @brief Calculate the @a binned distance and combined weight between this and four other CompactCoordinatesXYZW.
             *        This leverages more efficient SIMD instructions by using a 128-bit registers (SSE).
             */
            xyzw::QuadEvaluatedResultRounded evaluate_rounded(const CompactCoordinatesXYZW& v1, const CompactCoordinatesXYZW& v2, const CompactCoordinatesXYZW& v3, const CompactCoordinatesXYZW& v4) const noexcept;

            /**
             * @brief Calculate the distance and combined weight between this and four other CompactCoordinatesXYZW.
             *        This leverages more efficient SIMD instructions by using a 128-bit registers (SSE).
             */
            xyzw::QuadEvaluatedResult evaluate(const CompactCoordinatesXYZW& v1, const CompactCoordinatesXYZW& v2, const CompactCoordinatesXYZW& v3, const CompactCoordinatesXYZW& v4) const noexcept;

            /**
             * @brief Calculate the @a binned distance and combined weight between this and four other CompactCoordinatesXYZW.
             *        This leverages more efficient SIMD instructions by using either two 128-bit registers (SSE) or one 256-bit register (AVX).
             */
            xyzw::OctoEvaluatedResultRounded evaluate_rounded(
                const CompactCoordinatesXYZW& v1, const CompactCoordinatesXYZW& v2, const CompactCoordinatesXYZW& v3, const CompactCoordinatesXYZW& v4,
                const CompactCoordinatesXYZW& v5, const CompactCoordinatesXYZW& v6, const CompactCoordinatesXYZW& v7, const CompactCoordinatesXYZW& v8
            ) const noexcept;

            /**
             * @brief Calculate the distance and combined weight between this and four other CompactCoordinatesXYZW.
             *        This leverages more efficient SIMD instructions by using either two 128-bit registers (SSE) or one 256-bit register (AVX).
             */
            xyzw::OctoEvaluatedResult evaluate(
                const CompactCoordinatesXYZW& v1, const CompactCoordinatesXYZW& v2, const CompactCoordinatesXYZW& v3, const CompactCoordinatesXYZW& v4,
                const CompactCoordinatesXYZW& v5, const CompactCoordinatesXYZW& v6, const CompactCoordinatesXYZW& v7, const CompactCoordinatesXYZW& v8
            ) const noexcept;

            union {
                struct {Vector3<float> pos; float w;} value;
                std::array<float, 4> data;
            };

        protected:
            xyzw::EvaluatedResultRounded evaluate_rounded_scalar(const CompactCoordinatesXYZW& other) const noexcept;
            xyzw::EvaluatedResult evaluate_scalar(const CompactCoordinatesXYZW& other) const noexcept;

            xyzw::QuadEvaluatedResultRounded evaluate_rounded_scalar(const CompactCoordinatesXYZW& v1, const CompactCoordinatesXYZW& v2, const CompactCoordinatesXYZW& v3, const CompactCoordinatesXYZW& v4) const noexcept;
            xyzw::QuadEvaluatedResult evaluate_scalar(const CompactCoordinatesXYZW& v1, const CompactCoordinatesXYZW& v2, const CompactCoordinatesXYZW& v3, const CompactCoordinatesXYZW& v4) const noexcept;

            xyzw::OctoEvaluatedResultRounded evaluate_rounded_scalar(
                const CompactCoordinatesXYZW& v1, const CompactCoordinatesXYZW& v2, const CompactCoordinatesXYZW& v3, const CompactCoordinatesXYZW& v4,
                const CompactCoordinatesXYZW& v5, const CompactCoordinatesXYZW& v6, const CompactCoordinatesXYZW& v7, const CompactCoordinatesXYZW& v8
            ) const noexcept;
            xyzw::OctoEvaluatedResult evaluate_scalar(
                const CompactCoordinatesXYZW& v1, const CompactCoordinatesXYZW& v2, const CompactCoordinatesXYZW& v3, const CompactCoordinatesXYZW& v4,
                const CompactCoordinatesXYZW& v5, const CompactCoordinatesXYZW& v6, const CompactCoordinatesXYZW& v7, const CompactCoordinatesXYZW& v8
            ) const noexcept;

            #if defined __SSE2__
                xyzw::EvaluatedResultRounded evaluate_rounded_sse(const CompactCoordinatesXYZW& other) const noexcept;
                xyzw::EvaluatedResult evaluate_sse(const CompactCoordinatesXYZW& other) const noexcept;

                xyzw::QuadEvaluatedResultRounded evaluate_rounded_sse(const CompactCoordinatesXYZW& v1, const CompactCoordinatesXYZW& v2, const CompactCoordinatesXYZW& v3, const CompactCoordinatesXYZW& v4) const noexcept;
                xyzw::QuadEvaluatedResult evaluate_sse(const CompactCoordinatesXYZW& v1, const CompactCoordinatesXYZW& v2, const CompactCoordinatesXYZW& v3, const CompactCoordinatesXYZW& v4) const noexcept;

                xyzw::OctoEvaluatedResultRounded evaluate_rounded_sse(
                    const CompactCoordinatesXYZW& v1, const CompactCoordinatesXYZW& v2, const CompactCoordinatesXYZW& v3, const CompactCoordinatesXYZW& v4,
                    const CompactCoordinatesXYZW& v5, const CompactCoordinatesXYZW& v6, const CompactCoordinatesXYZW& v7, const CompactCoordinatesXYZW& v8
                ) const noexcept;
                xyzw::OctoEvaluatedResult evaluate_sse(
                    const CompactCoordinatesXYZW& v1, const CompactCoordinatesXYZW& v2, const CompactCoordinatesXYZW& v3, const CompactCoordinatesXYZW& v4,
                    const CompactCoordinatesXYZW& v5, const CompactCoordinatesXYZW& v6, const CompactCoordinatesXYZW& v7, const CompactCoordinatesXYZW& v8
                ) const noexcept;
            #endif

            #if defined __AVX__
                xyzw::EvaluatedResultRounded evaluate_rounded_avx(const CompactCoordinatesXYZW& other) const noexcept;
                xyzw::EvaluatedResult evaluate_avx(const CompactCoordinatesXYZW& other) const noexcept;

                xyzw::QuadEvaluatedResultRounded evaluate_rounded_avx(const CompactCoordinatesXYZW& v1, const CompactCoordinatesXYZW& v2, const CompactCoordinatesXYZW& v3, const CompactCoordinatesXYZW& v4) const noexcept;
                xyzw::QuadEvaluatedResult evaluate_avx(const CompactCoordinatesXYZW& v1, const CompactCoordinatesXYZW& v2, const CompactCoordinatesXYZW& v3, const CompactCoordinatesXYZW& v4) const noexcept;

                xyzw::OctoEvaluatedResultRounded evaluate_rounded_avx(
                    const CompactCoordinatesXYZW& v1, const CompactCoordinatesXYZW& v2, const CompactCoordinatesXYZW& v3, const CompactCoordinatesXYZW& v4,
                    const CompactCoordinatesXYZW& v5, const CompactCoordinatesXYZW& v6, const CompactCoordinatesXYZW& v7, const CompactCoordinatesXYZW& v8
                ) const noexcept;
                xyzw::OctoEvaluatedResult evaluate_avx(
                    const CompactCoordinatesXYZW& v1, const CompactCoordinatesXYZW& v2, const CompactCoordinatesXYZW& v3, const CompactCoordinatesXYZW& v4,
                    const CompactCoordinatesXYZW& v5, const CompactCoordinatesXYZW& v6, const CompactCoordinatesXYZW& v7, const CompactCoordinatesXYZW& v8
                ) const noexcept;
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

// undefine SSE2 & AVX for MacOS
#if defined __APPLE__
    #undef __AVX__
    #undef __SSE2__
    #undef __SSE4_2__
#endif

template<bool vbw>
inline ausaxs::hist::detail::xyzw::EvaluatedResult ausaxs::hist::detail::CompactCoordinatesXYZW<vbw>::evaluate(const CompactCoordinatesXYZW& other) const noexcept {
    #if defined __SSE2__
        return evaluate_sse(other);
    #else
        return evaluate_scalar(other);
    #endif
}

template<bool vbw>
inline ausaxs::hist::detail::xyzw::EvaluatedResultRounded ausaxs::hist::detail::CompactCoordinatesXYZW<vbw>::evaluate_rounded(const CompactCoordinatesXYZW& other) const noexcept {
    #if defined __SSE2__
        return evaluate_rounded_sse(other);
    #else
        return evaluate_rounded_scalar(other);
    #endif
}

template<bool vbw>
inline ausaxs::hist::detail::xyzw::QuadEvaluatedResult ausaxs::hist::detail::CompactCoordinatesXYZW<vbw>::evaluate(
    const CompactCoordinatesXYZW& v1, const CompactCoordinatesXYZW& v2, const CompactCoordinatesXYZW& v3, const CompactCoordinatesXYZW& v4
) const noexcept {
    #if defined __AVX__
        return evaluate_avx(v1, v2, v3, v4);
    #elif defined __SSE2__
        return evaluate_sse(v1, v2, v3, v4);
    #else
        return evaluate_scalar(v1, v2, v3, v4);
    #endif
}

template<bool vbw>
inline ausaxs::hist::detail::xyzw::QuadEvaluatedResultRounded ausaxs::hist::detail::CompactCoordinatesXYZW<vbw>::evaluate_rounded(
    const CompactCoordinatesXYZW& v1, const CompactCoordinatesXYZW& v2, const CompactCoordinatesXYZW& v3, const CompactCoordinatesXYZW& v4
) const noexcept {
    #if defined __AVX__
        return evaluate_rounded_avx(v1, v2, v3, v4);
    #elif defined __SSE2__
        return evaluate_rounded_sse(v1, v2, v3, v4);
    #else
        return evaluate_rounded_scalar(v1, v2, v3, v4);
    #endif
}

template<bool vbw>
inline ausaxs::hist::detail::xyzw::OctoEvaluatedResult ausaxs::hist::detail::CompactCoordinatesXYZW<vbw>::evaluate(
    const CompactCoordinatesXYZW& v1, const CompactCoordinatesXYZW& v2, const CompactCoordinatesXYZW& v3, const CompactCoordinatesXYZW& v4, 
    const CompactCoordinatesXYZW& v5, const CompactCoordinatesXYZW& v6, const CompactCoordinatesXYZW& v7, const CompactCoordinatesXYZW& v8
) const noexcept {
    #if defined __AVX__
        return evaluate_avx(v1, v2, v3, v4, v5, v6, v7, v8);
    #elif defined __SSE2__
        return evaluate_sse(v1, v2, v3, v4, v5, v6, v7, v8);
    #else
        return evaluate_scalar(v1, v2, v3, v4, v5, v6, v7, v8);
    #endif
}

template<bool vbw>
inline ausaxs::hist::detail::xyzw::OctoEvaluatedResultRounded ausaxs::hist::detail::CompactCoordinatesXYZW<vbw>::evaluate_rounded(
    const CompactCoordinatesXYZW& v1, const CompactCoordinatesXYZW& v2, const CompactCoordinatesXYZW& v3, const CompactCoordinatesXYZW& v4, 
    const CompactCoordinatesXYZW& v5, const CompactCoordinatesXYZW& v6, const CompactCoordinatesXYZW& v7, const CompactCoordinatesXYZW& v8
) const noexcept {
    #if defined __AVX__
        return evaluate_rounded_avx(v1, v2, v3, v4, v5, v6, v7, v8);
    #elif defined __SSE2__
        return evaluate_rounded_sse(v1, v2, v3, v4, v5, v6, v7, v8);
    #else
        return evaluate_rounded_scalar(v1, v2, v3, v4, v5, v6, v7, v8);
    #endif
}

template<bool vbw>
inline ausaxs::hist::detail::xyzw::EvaluatedResult ausaxs::hist::detail::CompactCoordinatesXYZW<vbw>::evaluate_scalar(const CompactCoordinatesXYZW& other) const noexcept {
    float dist = std::sqrt(squared_dot_product(this->data.data(), other.data.data()));
    return xyzw::EvaluatedResult(dist, value.w*other.value.w);
}

template<bool vbw>
inline ausaxs::hist::detail::xyzw::EvaluatedResultRounded ausaxs::hist::detail::CompactCoordinatesXYZW<vbw>::evaluate_rounded_scalar(
    const CompactCoordinatesXYZW& other
) const noexcept {
    int32_t dist = std::round(get_inv_width()*std::sqrt(squared_dot_product(this->data.data(), other.data.data())));
    return xyzw::EvaluatedResultRounded(dist, value.w*other.value.w);
}

template<bool vbw>
inline ausaxs::hist::detail::xyzw::QuadEvaluatedResult ausaxs::hist::detail::CompactCoordinatesXYZW<vbw>::evaluate_scalar(
    const CompactCoordinatesXYZW& v1, const CompactCoordinatesXYZW& v2, const CompactCoordinatesXYZW& v3, const CompactCoordinatesXYZW& v4
) const noexcept {
    float dx1 = std::sqrt(squared_dot_product(this->data.data(), v1.data.data()));
    float dx2 = std::sqrt(squared_dot_product(this->data.data(), v2.data.data()));
    float dx3 = std::sqrt(squared_dot_product(this->data.data(), v3.data.data()));
    float dx4 = std::sqrt(squared_dot_product(this->data.data(), v4.data.data()));
    return xyzw::QuadEvaluatedResult(
        std::array<float, 4>{dx1, dx2, dx3, dx4},
        std::array<float, 4>{value.w*v1.value.w, value.w*v2.value.w, value.w*v3.value.w, value.w*v4.value.w}
    );
}

template<bool vbw>
inline ausaxs::hist::detail::xyzw::QuadEvaluatedResultRounded ausaxs::hist::detail::CompactCoordinatesXYZW<vbw>::evaluate_rounded_scalar(
    const CompactCoordinatesXYZW& v1, const CompactCoordinatesXYZW& v2, const CompactCoordinatesXYZW& v3, const CompactCoordinatesXYZW& v4
) const noexcept {
    int32_t dx1 = std::round(get_inv_width()*std::sqrt(squared_dot_product(this->data.data(), v1.data.data())));
    int32_t dx2 = std::round(get_inv_width()*std::sqrt(squared_dot_product(this->data.data(), v2.data.data())));
    int32_t dx3 = std::round(get_inv_width()*std::sqrt(squared_dot_product(this->data.data(), v3.data.data())));
    int32_t dx4 = std::round(get_inv_width()*std::sqrt(squared_dot_product(this->data.data(), v4.data.data())));
    return xyzw::QuadEvaluatedResultRounded(
        std::array<int32_t, 4>{dx1, dx2, dx3, dx4},
        std::array<float, 4>{value.w*v1.value.w, value.w*v2.value.w, value.w*v3.value.w, value.w*v4.value.w}
    );
}

template<bool vbw>
inline ausaxs::hist::detail::xyzw::OctoEvaluatedResult ausaxs::hist::detail::CompactCoordinatesXYZW<vbw>::evaluate_scalar(
    const CompactCoordinatesXYZW& v1, const CompactCoordinatesXYZW& v2, const CompactCoordinatesXYZW& v3, const CompactCoordinatesXYZW& v4, 
    const CompactCoordinatesXYZW& v5, const CompactCoordinatesXYZW& v6, const CompactCoordinatesXYZW& v7, const CompactCoordinatesXYZW& v8
) const noexcept {
    float dx1 = std::sqrt(squared_dot_product(this->data.data(), v1.data.data()));
    float dx2 = std::sqrt(squared_dot_product(this->data.data(), v2.data.data()));
    float dx3 = std::sqrt(squared_dot_product(this->data.data(), v3.data.data()));
    float dx4 = std::sqrt(squared_dot_product(this->data.data(), v4.data.data()));
    float dx5 = std::sqrt(squared_dot_product(this->data.data(), v5.data.data()));
    float dx6 = std::sqrt(squared_dot_product(this->data.data(), v6.data.data()));
    float dx7 = std::sqrt(squared_dot_product(this->data.data(), v7.data.data()));
    float dx8 = std::sqrt(squared_dot_product(this->data.data(), v8.data.data()));
    return xyzw::OctoEvaluatedResult(
        std::array<float, 8>{dx1, dx2, dx3, dx4, dx5, dx6, dx7, dx8},
        std::array<float, 8>{value.w*v1.value.w, value.w*v2.value.w, value.w*v3.value.w, value.w*v4.value.w, value.w*v5.value.w, value.w*v6.value.w, value.w*v7.value.w, value.w*v8.value.w}
    );
}

template<bool vbw>
inline ausaxs::hist::detail::xyzw::OctoEvaluatedResultRounded ausaxs::hist::detail::CompactCoordinatesXYZW<vbw>::evaluate_rounded_scalar(
    const CompactCoordinatesXYZW& v1, const CompactCoordinatesXYZW& v2, const CompactCoordinatesXYZW& v3, const CompactCoordinatesXYZW& v4, 
    const CompactCoordinatesXYZW& v5, const CompactCoordinatesXYZW& v6, const CompactCoordinatesXYZW& v7, const CompactCoordinatesXYZW& v8
) const noexcept {
    int32_t dx1 = std::round(get_inv_width()*std::sqrt(squared_dot_product(this->data.data(), v1.data.data())));
    int32_t dx2 = std::round(get_inv_width()*std::sqrt(squared_dot_product(this->data.data(), v2.data.data())));
    int32_t dx3 = std::round(get_inv_width()*std::sqrt(squared_dot_product(this->data.data(), v3.data.data())));
    int32_t dx4 = std::round(get_inv_width()*std::sqrt(squared_dot_product(this->data.data(), v4.data.data())));
    int32_t dx5 = std::round(get_inv_width()*std::sqrt(squared_dot_product(this->data.data(), v5.data.data())));
    int32_t dx6 = std::round(get_inv_width()*std::sqrt(squared_dot_product(this->data.data(), v6.data.data())));
    int32_t dx7 = std::round(get_inv_width()*std::sqrt(squared_dot_product(this->data.data(), v7.data.data())));
    int32_t dx8 = std::round(get_inv_width()*std::sqrt(squared_dot_product(this->data.data(), v8.data.data())));
    return xyzw::OctoEvaluatedResultRounded(
        std::array<int32_t, 8>{dx1, dx2, dx3, dx4, dx5, dx6, dx7, dx8},
        std::array<float, 8>{value.w*v1.value.w, value.w*v2.value.w, value.w*v3.value.w, value.w*v4.value.w, value.w*v5.value.w, value.w*v6.value.w, value.w*v7.value.w, value.w*v8.value.w}
    );
}

#if defined __SSE2__
    template<bool vbw>
    inline ausaxs::hist::detail::xyzw::EvaluatedResult ausaxs::hist::detail::CompactCoordinatesXYZW<vbw>::evaluate_sse(
        const CompactCoordinatesXYZW& other
    ) const noexcept {
        __m128 dist2 = squared_dot_product(this->data.data(), other.data.data(), OutputControl::ALL);
        __m128 dist_sqrt = _mm_sqrt_ps(dist2);
        float dist = _mm_cvtss_f32(dist_sqrt);
        return xyzw::EvaluatedResult(dist, this->value.w*other.value.w);
    }

    template<bool vbw>
    inline ausaxs::hist::detail::xyzw::EvaluatedResultRounded ausaxs::hist::detail::CompactCoordinatesXYZW<vbw>::evaluate_rounded_sse(
        const CompactCoordinatesXYZW& other
    ) const noexcept {
        __m128 dist2 = squared_dot_product(this->data.data(), other.data.data(), OutputControl::ALL);
        __m128 dist_sqrt = _mm_sqrt_ps(dist2);
        int32_t dist_bin = std::round(get_inv_width()*_mm_cvtss_f32(dist_sqrt));
        return xyzw::EvaluatedResultRounded(dist_bin, this->value.w*other.value.w);
    }

    template<bool vbw>
    inline ausaxs::hist::detail::xyzw::QuadEvaluatedResult ausaxs::hist::detail::CompactCoordinatesXYZW<vbw>::evaluate_sse(
        const CompactCoordinatesXYZW& v1, const CompactCoordinatesXYZW& v2, const CompactCoordinatesXYZW& v3, const CompactCoordinatesXYZW& v4
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

        __m128 w1 = _mm_set_ps1(value.w);
        __m128 w2 = _mm_set_ps(v4.value.w, v3.value.w, v2.value.w, v1.value.w);
        __m128 weights = _mm_mul_ps(w1, w2);

        xyzw::QuadEvaluatedResult result;
        _mm_store_ps(reinterpret_cast<float*>(result.distances.data()), dist_sqrt);         // efficient store of distances
        _mm_store_ps(reinterpret_cast<float*>(result.weights.data()), weights);             // efficient store of weights
        return result;
    }

    template<bool vbw>
    inline ausaxs::hist::detail::xyzw::QuadEvaluatedResultRounded ausaxs::hist::detail::CompactCoordinatesXYZW<vbw>::evaluate_rounded_sse(
        const CompactCoordinatesXYZW& v1, const CompactCoordinatesXYZW& v2, const CompactCoordinatesXYZW& v3, const CompactCoordinatesXYZW& v4
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
        __m128i dist_bin = _mm_cvtps_epi32(dist_binf);                          // cast to int

        __m128 w1 = _mm_set_ps1(value.w);
        __m128 w2 = _mm_set_ps(v4.value.w, v3.value.w, v2.value.w, v1.value.w);
        __m128 weights = _mm_mul_ps(w1, w2);

        xyzw::QuadEvaluatedResultRounded result;
        _mm_store_si128(reinterpret_cast<__m128i*>(result.distances.data()), dist_bin); // efficient store of bins
        _mm_store_ps(reinterpret_cast<float*>(result.weights.data()), weights);         // efficient store of weights
        return result;
    }

    template<bool vbw>
    inline ausaxs::hist::detail::xyzw::OctoEvaluatedResult ausaxs::hist::detail::CompactCoordinatesXYZW<vbw>::evaluate_sse(
        const CompactCoordinatesXYZW& v1, const CompactCoordinatesXYZW& v2, const CompactCoordinatesXYZW& v3, const CompactCoordinatesXYZW& v4, 
        const CompactCoordinatesXYZW& v5, const CompactCoordinatesXYZW& v6, const CompactCoordinatesXYZW& v7, const CompactCoordinatesXYZW& v8
    ) const noexcept {
        xyzw::OctoEvaluatedResult result;
        {   // first four
            __m128 dist2_1 = squared_dot_product(this->data.data(), v1.data.data(), OutputControl::FIRST);
            __m128 dist2_2 = squared_dot_product(this->data.data(), v2.data.data(), OutputControl::SECOND);
            __m128 dist2_3 = squared_dot_product(this->data.data(), v3.data.data(), OutputControl::THIRD);
            __m128 dist2_4 = squared_dot_product(this->data.data(), v4.data.data(), OutputControl::FOURTH);

            __m128 dist2 = _mm_add_ps(dist2_1, dist2_2);
            dist2 = _mm_add_ps(dist2, dist2_3);
            dist2 = _mm_add_ps(dist2, dist2_4);
            __m128 dist_sqrt = _mm_sqrt_ps(dist2);

            __m128 w1 = _mm_set_ps1(value.w);
            __m128 w2 = _mm_set_ps(v4.value.w, v3.value.w, v2.value.w, v1.value.w);
            __m128 weights = _mm_mul_ps(w1, w2);

            _mm_store_ps(reinterpret_cast<float*>(result.distances.data()), dist_sqrt);
            _mm_store_ps(reinterpret_cast<float*>(result.weights.data()), weights);
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

            __m128 w1 = _mm_set_ps1(value.w);
            __m128 w2 = _mm_set_ps(v8.value.w, v7.value.w, v6.value.w, v5.value.w);
            __m128 weights = _mm_mul_ps(w1, w2);

            _mm_store_ps(reinterpret_cast<float*>(result.distances.data()+4), dist_sqrt);
            _mm_store_ps(reinterpret_cast<float*>(result.weights.data()+4), weights);
        }
        return result;
    }

    template<bool vbw>
    inline ausaxs::hist::detail::xyzw::OctoEvaluatedResultRounded ausaxs::hist::detail::CompactCoordinatesXYZW<vbw>::evaluate_rounded_sse(
        const CompactCoordinatesXYZW& v1, const CompactCoordinatesXYZW& v2, const CompactCoordinatesXYZW& v3, const CompactCoordinatesXYZW& v4, 
        const CompactCoordinatesXYZW& v5, const CompactCoordinatesXYZW& v6, const CompactCoordinatesXYZW& v7, const CompactCoordinatesXYZW& v8
    ) const noexcept {
        xyzw::OctoEvaluatedResultRounded result;
        {   // first four
            __m128 dist2_1 = squared_dot_product(this->data.data(), v1.data.data(), OutputControl::FIRST);
            __m128 dist2_2 = squared_dot_product(this->data.data(), v2.data.data(), OutputControl::SECOND);
            __m128 dist2_3 = squared_dot_product(this->data.data(), v3.data.data(), OutputControl::THIRD);
            __m128 dist2_4 = squared_dot_product(this->data.data(), v4.data.data(), OutputControl::FOURTH);

            __m128 dist2 = _mm_add_ps(dist2_1, dist2_2);
            dist2 = _mm_add_ps(dist2, dist2_3);
            dist2 = _mm_add_ps(dist2, dist2_4);
            __m128 dist_sqrt = _mm_sqrt_ps(dist2);
            __m128 dist_binf = _mm_mul_ps(dist_sqrt, _mm_set_ps1(get_inv_width())); // multiply by the inverse width
            __m128i dist_bin = _mm_cvtps_epi32(dist_binf);                          // cast to int

            __m128 w1 = _mm_set_ps1(value.w);
            __m128 w2 = _mm_set_ps(v4.value.w, v3.value.w, v2.value.w, v1.value.w);
            __m128 weights = _mm_mul_ps(w1, w2);

            _mm_store_si128(reinterpret_cast<__m128i*>(result.distances.data()), dist_bin);
            _mm_store_ps(reinterpret_cast<float*>(result.weights.data()), weights);
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
            __m128 dist_binf = _mm_mul_ps(dist_sqrt, _mm_set_ps1(get_inv_width())); // multiply by the inverse width
            __m128i dist_bin = _mm_cvtps_epi32(dist_binf);                          // cast to int

            __m128 w1 = _mm_set_ps1(value.w);
            __m128 w2 = _mm_set_ps(v8.value.w, v7.value.w, v6.value.w, v5.value.w);
            __m128 weights = _mm_mul_ps(w1, w2);

            _mm_store_si128(reinterpret_cast<__m128i*>(result.distances.data()+4), dist_bin);
            _mm_store_ps(reinterpret_cast<float*>(result.weights.data()+4), weights);
        }
        return result;
    }
#endif

#if defined __AVX__
    template<bool vbw>
    inline ausaxs::hist::detail::xyzw::EvaluatedResult ausaxs::hist::detail::CompactCoordinatesXYZW<vbw>::evaluate_avx(
        const CompactCoordinatesXYZW& other
    ) const noexcept {
        return evaluate_sse(other); // no way to optimize a single evaluation with AVX
    }

    template<bool vbw>
    inline ausaxs::hist::detail::xyzw::EvaluatedResultRounded ausaxs::hist::detail::CompactCoordinatesXYZW<vbw>::evaluate_rounded_avx(
        const CompactCoordinatesXYZW& other
    ) const noexcept {
        return evaluate_rounded_sse(other); // no way to optimize a single evaluation with AVX
    }

    template<bool vbw>
    inline ausaxs::hist::detail::xyzw::QuadEvaluatedResult ausaxs::hist::detail::CompactCoordinatesXYZW<vbw>::evaluate_avx(
        const CompactCoordinatesXYZW& v1, const CompactCoordinatesXYZW& v2, const CompactCoordinatesXYZW& v3, const CompactCoordinatesXYZW& v4
    ) const noexcept {
        __m256 dist2_1 = squared_dot_product(this->data.data(), v1.data.data(), v3.data.data(), OutputControl::FIRST); // |Δx1^2|0    |0    |0    |Δx3^2|0    |0    |0    |
        __m256 dist2_2 = squared_dot_product(this->data.data(), v2.data.data(), v4.data.data(), OutputControl::SECOND);// |0    |Δx2^2|0    |0    |0    |Δx4^2|0    |0    |

        __m256 dist2_256 = _mm256_add_ps(dist2_1, dist2_2);                     // |Δx1^2|Δx2^2|0    |0    |Δx3^2|Δx4^2|0    |0    |
        __m256 dist2_256_shuffled = _mm256_permute_ps(dist2_256, 0b01001110);   // |0    |0    |Δx3^2|Δx4^2|0    |0    |Δx1^2|Δx2^2|
        __m128 dist2_128_lower = _mm256_extractf128_ps(dist2_256, 0);           // |Δx1^2|Δx2^2|0    |0    |
        __m128 dist2_128_upper = _mm256_extractf128_ps(dist2_256_shuffled, 1);  // |0    |0    |Δx3^2|Δx4^2|
        __m128 dist2 = _mm_add_ps(dist2_128_lower, dist2_128_upper);            // |Δx1^2|Δx2^2|Δx3^2|Δx4^2|
        __m128 dist_sqrt = _mm_sqrt_ps(dist2);

        __m128 w1 = _mm_set1_ps(value.w);
        __m128 w2 = _mm_set_ps(v4.value.w, v3.value.w, v2.value.w, v1.value.w);
        __m128 weights = _mm_mul_ps(w1, w2);

        xyzw::QuadEvaluatedResult result;
        _mm_store_ps(reinterpret_cast<float*>(result.distances.data()), dist_sqrt);         // efficient store of distances
        _mm_store_ps(reinterpret_cast<float*>(result.weights.data()), weights);             // efficient store of weights
        return result;
    }

    template<bool vbw>
    inline ausaxs::hist::detail::xyzw::QuadEvaluatedResultRounded ausaxs::hist::detail::CompactCoordinatesXYZW<vbw>::evaluate_rounded_avx(
        const CompactCoordinatesXYZW& v1, const CompactCoordinatesXYZW& v2, const CompactCoordinatesXYZW& v3, const CompactCoordinatesXYZW& v4
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
        __m128i dist_bin = _mm_cvtps_epi32(dist_binf);                          // cast to int

        __m128 w1 = _mm_set1_ps(value.w);
        __m128 w2 = _mm_set_ps(v4.value.w, v3.value.w, v2.value.w, v1.value.w);
        __m128 weights = _mm_mul_ps(w1, w2);

        xyzw::QuadEvaluatedResultRounded result;
        _mm_store_si128(reinterpret_cast<__m128i*>(result.distances.data()), dist_bin);  // efficient store of bins
        _mm_store_ps(reinterpret_cast<float*>(result.weights.data()), weights);          // efficient store of weights
        return result;
    }

    template<bool vbw>
    inline ausaxs::hist::detail::xyzw::OctoEvaluatedResult ausaxs::hist::detail::CompactCoordinatesXYZW<vbw>::evaluate_avx(
        const CompactCoordinatesXYZW& v1, const CompactCoordinatesXYZW& v2, const CompactCoordinatesXYZW& v3, const CompactCoordinatesXYZW& v4, 
        const CompactCoordinatesXYZW& v5, const CompactCoordinatesXYZW& v6, const CompactCoordinatesXYZW& v7, const CompactCoordinatesXYZW& v8
    ) const noexcept {
        __m256 dist2_1 = squared_dot_product(this->data.data(), v1.data.data(), v5.data.data(), OutputControl::FIRST); // |Δx1^2|0    |0    |0    |Δx5^2|0    |0    |0    |
        __m256 dist2_2 = squared_dot_product(this->data.data(), v2.data.data(), v6.data.data(), OutputControl::SECOND);// |0    |Δx2^2|0    |0    |0    |Δx6^2|0    |0    |
        __m256 dist2_3 = squared_dot_product(this->data.data(), v3.data.data(), v7.data.data(), OutputControl::THIRD); // |0    |0    |Δx3^2|0    |0    |0    |Δx7^2|0    |
        __m256 dist2_4 = squared_dot_product(this->data.data(), v4.data.data(), v8.data.data(), OutputControl::FOURTH);// |0    |0    |0    |Δx4^2|0    |0    |0    |Δx8^2|

        __m256 dist2_256 = _mm256_add_ps(dist2_1, dist2_2);                                                            // |Δx1^2|Δx2^2|0    |0    |Δx5^2|Δx6^2|0    |0    |
        dist2_256 = _mm256_add_ps(dist2_256, dist2_3);                                                                 // |Δx1^2|Δx2^2|Δx3^2|0    |Δx5^2|Δx6^2|Δx7^2|0    |
        dist2_256 = _mm256_add_ps(dist2_256, dist2_4);                                                                 // |Δx1^2|Δx2^2|Δx3^2|Δx4^2|Δx5^2|Δx6^2|Δx7^2|Δx8^2|
        __m256 dist_sqrt = _mm256_sqrt_ps(dist2_256);

        __m256 w1 = _mm256_set1_ps(value.w);
        __m256 w2 = _mm256_set_ps(v8.value.w, v7.value.w, v6.value.w, v5.value.w, v4.value.w, v3.value.w, v2.value.w, v1.value.w);
        __m256 weights = _mm256_mul_ps(w1, w2);

        xyzw::OctoEvaluatedResult result;
        _mm256_store_ps(reinterpret_cast<float*>(result.distances.data()), dist_sqrt);          // efficient store of distances
        _mm256_store_ps(reinterpret_cast<float*>(result.weights.data()), weights);              // efficient store of weights
        return result;
    }

    template<bool vbw>
    inline ausaxs::hist::detail::xyzw::OctoEvaluatedResultRounded ausaxs::hist::detail::CompactCoordinatesXYZW<vbw>::evaluate_rounded_avx(
        const CompactCoordinatesXYZW& v1, const CompactCoordinatesXYZW& v2, const CompactCoordinatesXYZW& v3, const CompactCoordinatesXYZW& v4, 
        const CompactCoordinatesXYZW& v5, const CompactCoordinatesXYZW& v6, const CompactCoordinatesXYZW& v7, const CompactCoordinatesXYZW& v8
    ) const noexcept {
        __m256 dist2_1 = squared_dot_product(this->data.data(), v1.data.data(), v5.data.data(), OutputControl::FIRST); // |Δx1^2|0    |0    |0    |Δx5^2|0    |0    |0    |
        __m256 dist2_2 = squared_dot_product(this->data.data(), v2.data.data(), v6.data.data(), OutputControl::SECOND);// |0    |Δx2^2|0    |0    |0    |Δx6^2|0    |0    |
        __m256 dist2_3 = squared_dot_product(this->data.data(), v3.data.data(), v7.data.data(), OutputControl::THIRD); // |0    |0    |Δx3^2|0    |0    |0    |Δx7^2|0    |
        __m256 dist2_4 = squared_dot_product(this->data.data(), v4.data.data(), v8.data.data(), OutputControl::FOURTH);// |0    |0    |0    |Δx4^2|0    |0    |0    |Δx8^2|

        __m256 dist2_256 = _mm256_add_ps(dist2_1, dist2_2);                                                            // |Δx1^2|Δx2^2|0    |0    |Δx5^2|Δx6^2|0    |0    |
        dist2_256 = _mm256_add_ps(dist2_256, dist2_3);                                                                 // |Δx1^2|Δx2^2|Δx3^2|0    |Δx5^2|Δx6^2|Δx7^2|0    |
        dist2_256 = _mm256_add_ps(dist2_256, dist2_4);                                                                 // |Δx1^2|Δx2^2|Δx3^2|Δx4^2|Δx5^2|Δx6^2|Δx7^2|Δx8^2|
        __m256 dist_sqrt = _mm256_sqrt_ps(dist2_256);
        __m256 dist_binf = _mm256_mul_ps(dist_sqrt, _mm256_set1_ps(get_inv_width()));   // multiply by the inverse width
        __m256i dist_bin = _mm256_cvtps_epi32(dist_binf);                               // cast to int

        __m256 w1 = _mm256_set1_ps(value.w);
        __m256 w2 = _mm256_set_ps(v8.value.w, v7.value.w, v6.value.w, v5.value.w, v4.value.w, v3.value.w, v2.value.w, v1.value.w);
        __m256 weights = _mm256_mul_ps(w1, w2);

        xyzw::OctoEvaluatedResultRounded result;
        _mm256_store_si256(reinterpret_cast<__m256i*>(result.distances.data()), dist_bin);   // efficient store of bins
        _mm256_store_ps(reinterpret_cast<float*>(result.weights.data()), weights);           // efficient store of weights
        return result;
    }
#endif