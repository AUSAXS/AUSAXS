// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <math/Vector3.h>
#include <constants/Constants.h>
#include <settings/Flags.h>

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

namespace ausaxs::hist::detail {
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

    struct ConstantWidth {
        static consteval float get() {return 1./ausaxs::constants::axes::d_axis.width();}
    };

    struct VariableWidth {
        static float get() {return settings::flags::inv_bin_width;}
    };

    template<bool variable_bin_width>
    class CompactCoordinatesData {
        public:
            CompactCoordinatesData() noexcept = default;

            template<numeric T, numeric V>
            CompactCoordinatesData(const Vector3<T>& v, V w) noexcept : value{.pos={static_cast<float>(v.x()), static_cast<float>(v.y()), static_cast<float>(v.z())}, .w=static_cast<float>(w)} {}
            CompactCoordinatesData(const Vector3<float>& v, float w) noexcept : value{.pos=v, .w=w} {}

            constexpr static float get_inv_width() {
                if constexpr (variable_bin_width) {
                    return VariableWidth::get();
                } else {
                    return ConstantWidth::get();
                }
            }

            /**
             * @brief Calculate the @a binned distance and combined weight between this and a single other CompactCoordinatesData.
             */
            EvaluatedResultRounded evaluate_rounded(const CompactCoordinatesData& other) const noexcept;

            /**
             * @brief Calculate the distance and combined weight between this and a single other CompactCoordinatesData.
             */
            EvaluatedResult evaluate(const CompactCoordinatesData& other) const noexcept;

            /**
             * @brief Calculate the @a binned distance and combined weight between this and four other CompactCoordinatesData.
             *        This leverages more efficient SIMD instructions by using a 128-bit registers (SSE).
             */
            QuadEvaluatedResultRounded evaluate_rounded(const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4) const noexcept;

            /**
             * @brief Calculate the distance and combined weight between this and four other CompactCoordinatesData.
             *        This leverages more efficient SIMD instructions by using a 128-bit registers (SSE).
             */
            QuadEvaluatedResult evaluate(const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4) const noexcept;

            /**
             * @brief Calculate the @a binned distance and combined weight between this and four other CompactCoordinatesData.
             *        This leverages more efficient SIMD instructions by using either two 128-bit registers (SSE) or one 256-bit register (AVX).
             */
            OctoEvaluatedResultRounded evaluate_rounded(
                const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4,
                const CompactCoordinatesData& v5, const CompactCoordinatesData& v6, const CompactCoordinatesData& v7, const CompactCoordinatesData& v8
            ) const noexcept;

            /**
             * @brief Calculate the distance and combined weight between this and four other CompactCoordinatesData.
             *        This leverages more efficient SIMD instructions by using either two 128-bit registers (SSE) or one 256-bit register (AVX).
             */
            OctoEvaluatedResult evaluate(
                const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4,
                const CompactCoordinatesData& v5, const CompactCoordinatesData& v6, const CompactCoordinatesData& v7, const CompactCoordinatesData& v8
            ) const noexcept;

            union {
                struct {Vector3<float> pos; float w;} value;
                std::array<float, 4> data;
            };

        protected:
            EvaluatedResultRounded evaluate_rounded_scalar(const CompactCoordinatesData& other) const noexcept;
            EvaluatedResult evaluate_scalar(const CompactCoordinatesData& other) const noexcept;

            QuadEvaluatedResultRounded evaluate_rounded_scalar(const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4) const noexcept;
            QuadEvaluatedResult evaluate_scalar(const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4) const noexcept;

            OctoEvaluatedResultRounded evaluate_rounded_scalar(
                const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4,
                const CompactCoordinatesData& v5, const CompactCoordinatesData& v6, const CompactCoordinatesData& v7, const CompactCoordinatesData& v8
            ) const noexcept;
            OctoEvaluatedResult evaluate_scalar(
                const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4,
                const CompactCoordinatesData& v5, const CompactCoordinatesData& v6, const CompactCoordinatesData& v7, const CompactCoordinatesData& v8
            ) const noexcept;

            #if defined __SSE2__
                EvaluatedResultRounded evaluate_rounded_sse(const CompactCoordinatesData& other) const noexcept;
                EvaluatedResult evaluate_sse(const CompactCoordinatesData& other) const noexcept;

                QuadEvaluatedResultRounded evaluate_rounded_sse(const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4) const noexcept;
                QuadEvaluatedResult evaluate_sse(const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4) const noexcept;

                OctoEvaluatedResultRounded evaluate_rounded_sse(
                    const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4,
                    const CompactCoordinatesData& v5, const CompactCoordinatesData& v6, const CompactCoordinatesData& v7, const CompactCoordinatesData& v8
                ) const noexcept;
                OctoEvaluatedResult evaluate_sse(
                    const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4,
                    const CompactCoordinatesData& v5, const CompactCoordinatesData& v6, const CompactCoordinatesData& v7, const CompactCoordinatesData& v8
                ) const noexcept;
            #endif

            #if defined __AVX__
                EvaluatedResultRounded evaluate_rounded_avx(const CompactCoordinatesData& other) const noexcept;
                EvaluatedResult evaluate_avx(const CompactCoordinatesData& other) const noexcept;

                QuadEvaluatedResultRounded evaluate_rounded_avx(const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4) const noexcept;
                QuadEvaluatedResult evaluate_avx(const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4) const noexcept;

                OctoEvaluatedResultRounded evaluate_rounded_avx(
                    const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4,
                    const CompactCoordinatesData& v5, const CompactCoordinatesData& v6, const CompactCoordinatesData& v7, const CompactCoordinatesData& v8
                ) const noexcept;
                OctoEvaluatedResult evaluate_avx(
                    const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4,
                    const CompactCoordinatesData& v5, const CompactCoordinatesData& v6, const CompactCoordinatesData& v7, const CompactCoordinatesData& v8
                ) const noexcept;
            #endif
    };
    static_assert(sizeof(CompactCoordinatesData<true>) == 16,              "CompactCoordinatesData is not 16 bytes. This is required for aligning SIMD instructions.");
    static_assert(std::is_trivial_v<CompactCoordinatesData<true>>,         "CompactCoordinatesData is not trivial");
    static_assert(std::is_standard_layout_v<CompactCoordinatesData<true>>, "CompactCoordinatesData is not standard layout");
    static_assert(supports_nothrow_move_v<CompactCoordinatesData<true>>,   "CompactCoordinatesData should support nothrow move semantics.");
    static_assert(sizeof(CompactCoordinatesData<false>) == 16,             "CompactCoordinatesData is not 16 bytes. This is required for aligning SIMD instructions.");
    static_assert(std::is_trivial_v<CompactCoordinatesData<false>>,        "CompactCoordinatesData is not trivial");
    static_assert(std::is_standard_layout_v<CompactCoordinatesData<false>>,"CompactCoordinatesData is not standard layout");
    static_assert(supports_nothrow_move_v<CompactCoordinatesData<false>>,  "CompactCoordinatesData should support nothrow move semantics.");
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

namespace ausaxs::hist::detail {
    constexpr float inv_width = ausaxs::constants::axes::d_inv_width;
}

inline ausaxs::hist::detail::EvaluatedResult ausaxs::hist::detail::CompactCoordinatesData::evaluate(const CompactCoordinatesData& other) const  noexcept{
    #if defined __SSE2__
        return evaluate_sse(other);
    #else
        return evaluate_scalar(other);
    #endif
}

inline ausaxs::hist::detail::EvaluatedResultRounded ausaxs::hist::detail::CompactCoordinatesData::evaluate_rounded(const CompactCoordinatesData& other) const  noexcept{
    #if defined __SSE2__
        return evaluate_rounded_sse(other);
    #else
        return evaluate_rounded_scalar(other);
    #endif
}

template<bool vbw>
inline ausaxs::hist::detail::QuadEvaluatedResult ausaxs::hist::detail::CompactCoordinatesData<vbw>::evaluate(
    const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4
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
inline ausaxs::hist::detail::QuadEvaluatedResultRounded ausaxs::hist::detail::CompactCoordinatesData<vbw>::evaluate_rounded(
    const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4
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
inline ausaxs::hist::detail::OctoEvaluatedResult ausaxs::hist::detail::CompactCoordinatesData<vbw>::evaluate(
    const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4, 
    const CompactCoordinatesData& v5, const CompactCoordinatesData& v6, const CompactCoordinatesData& v7, const CompactCoordinatesData& v8
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
inline ausaxs::hist::detail::OctoEvaluatedResultRounded ausaxs::hist::detail::CompactCoordinatesData<vbw>::evaluate_rounded(
    const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4, 
    const CompactCoordinatesData& v5, const CompactCoordinatesData& v6, const CompactCoordinatesData& v7, const CompactCoordinatesData& v8
) const noexcept {
    #if defined __AVX__
        return evaluate_rounded_avx(v1, v2, v3, v4, v5, v6, v7, v8);
    #elif defined __SSE2__
        return evaluate_rounded_sse(v1, v2, v3, v4, v5, v6, v7, v8);
    #else
        return evaluate_rounded_scalar(v1, v2, v3, v4, v5, v6, v7, v8);
    #endif
}

static inline float squared_dot_product(const float* v1, const float* v2) noexcept {
    float dx = v1[0] - v2[0];
    float dy = v1[1] - v2[1];
    float dz = v1[2] - v2[2];
    return dx*dx + dy*dy + dz*dz;
}

inline ausaxs::hist::detail::EvaluatedResult ausaxs::hist::detail::CompactCoordinatesData::evaluate_scalar(
    const CompactCoordinatesData& other
) const noexcept {
    float dist = std::sqrt(squared_dot_product(this->data.data(), other.data.data()));
    return EvaluatedResult(dist, value.w*other.value.w);
}

inline ausaxs::hist::detail::EvaluatedResultRounded ausaxs::hist::detail::CompactCoordinatesData::evaluate_rounded_scalar(
    const CompactCoordinatesData& other
) const noexcept {
    int32_t dist = std::round(inv_width*std::sqrt(squared_dot_product(this->data.data(), other.data.data())));
    return EvaluatedResultRounded(dist, value.w*other.value.w);
}

template<bool vbw>
inline ausaxs::hist::detail::QuadEvaluatedResult ausaxs::hist::detail::CompactCoordinatesData<vbw>::evaluate_scalar(
    const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4
) const noexcept {
    float dx1 = std::sqrt(squared_dot_product(this->data.data(), v1.data.data()));
    float dx2 = std::sqrt(squared_dot_product(this->data.data(), v2.data.data()));
    float dx3 = std::sqrt(squared_dot_product(this->data.data(), v3.data.data()));
    float dx4 = std::sqrt(squared_dot_product(this->data.data(), v4.data.data()));
    return QuadEvaluatedResult(
        std::array<float, 4>{dx1, dx2, dx3, dx4},
        std::array<float, 4>{value.w*v1.value.w, value.w*v2.value.w, value.w*v3.value.w, value.w*v4.value.w}
    );
}

template<bool vbw>
inline ausaxs::hist::detail::QuadEvaluatedResultRounded ausaxs::hist::detail::CompactCoordinatesData<vbw>::evaluate_rounded_scalar(
    const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4
) const noexcept {
    int32_t dx1 = std::round(inv_width*std::sqrt(squared_dot_product(this->data.data(), v1.data.data())));
    int32_t dx2 = std::round(inv_width*std::sqrt(squared_dot_product(this->data.data(), v2.data.data())));
    int32_t dx3 = std::round(inv_width*std::sqrt(squared_dot_product(this->data.data(), v3.data.data())));
    int32_t dx4 = std::round(inv_width*std::sqrt(squared_dot_product(this->data.data(), v4.data.data())));
    return QuadEvaluatedResultRounded(
        std::array<int32_t, 4>{dx1, dx2, dx3, dx4},
        std::array<float, 4>{value.w*v1.value.w, value.w*v2.value.w, value.w*v3.value.w, value.w*v4.value.w}
    );
}

template<bool vbw>
inline ausaxs::hist::detail::OctoEvaluatedResult ausaxs::hist::detail::CompactCoordinatesData<vbw>::evaluate_scalar(
    const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4, 
    const CompactCoordinatesData& v5, const CompactCoordinatesData& v6, const CompactCoordinatesData& v7, const CompactCoordinatesData& v8
) const noexcept {
    float dx1 = std::sqrt(squared_dot_product(this->data.data(), v1.data.data()));
    float dx2 = std::sqrt(squared_dot_product(this->data.data(), v2.data.data()));
    float dx3 = std::sqrt(squared_dot_product(this->data.data(), v3.data.data()));
    float dx4 = std::sqrt(squared_dot_product(this->data.data(), v4.data.data()));
    float dx5 = std::sqrt(squared_dot_product(this->data.data(), v5.data.data()));
    float dx6 = std::sqrt(squared_dot_product(this->data.data(), v6.data.data()));
    float dx7 = std::sqrt(squared_dot_product(this->data.data(), v7.data.data()));
    float dx8 = std::sqrt(squared_dot_product(this->data.data(), v8.data.data()));
    return OctoEvaluatedResult(
        std::array<float, 8>{dx1, dx2, dx3, dx4, dx5, dx6, dx7, dx8},
        std::array<float, 8>{value.w*v1.value.w, value.w*v2.value.w, value.w*v3.value.w, value.w*v4.value.w, value.w*v5.value.w, value.w*v6.value.w, value.w*v7.value.w, value.w*v8.value.w}
    );
}

template<bool vbw>
inline ausaxs::hist::detail::OctoEvaluatedResultRounded ausaxs::hist::detail::CompactCoordinatesData<vbw>::evaluate_rounded_scalar(
    const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4, 
    const CompactCoordinatesData& v5, const CompactCoordinatesData& v6, const CompactCoordinatesData& v7, const CompactCoordinatesData& v8
) const noexcept {
    int32_t dx1 = std::round(inv_width*std::sqrt(squared_dot_product(this->data.data(), v1.data.data())));
    int32_t dx2 = std::round(inv_width*std::sqrt(squared_dot_product(this->data.data(), v2.data.data())));
    int32_t dx3 = std::round(inv_width*std::sqrt(squared_dot_product(this->data.data(), v3.data.data())));
    int32_t dx4 = std::round(inv_width*std::sqrt(squared_dot_product(this->data.data(), v4.data.data())));
    int32_t dx5 = std::round(inv_width*std::sqrt(squared_dot_product(this->data.data(), v5.data.data())));
    int32_t dx6 = std::round(inv_width*std::sqrt(squared_dot_product(this->data.data(), v6.data.data())));
    int32_t dx7 = std::round(inv_width*std::sqrt(squared_dot_product(this->data.data(), v7.data.data())));
    int32_t dx8 = std::round(inv_width*std::sqrt(squared_dot_product(this->data.data(), v8.data.data())));
    return OctoEvaluatedResultRounded(
        std::array<int32_t, 8>{dx1, dx2, dx3, dx4, dx5, dx6, dx7, dx8},
        std::array<float, 8>{value.w*v1.value.w, value.w*v2.value.w, value.w*v3.value.w, value.w*v4.value.w, value.w*v5.value.w, value.w*v6.value.w, value.w*v7.value.w, value.w*v8.value.w}
    );
}

#if defined __SSE2__
    #include <nmmintrin.h>
    namespace ausaxs::hist::detail {
        enum OutputControl : int8_t {
            ALL =    0b01111111,
            FIRST =  0b01110001,
            SECOND = 0b01110010,
            THIRD =  0b01110100,
            FOURTH = 0b01111000
        };

        /**
        * @brief Calculate the squared distance between two CompactCoordinatesData using 128-bit SSE2 instructions.
        */
        static inline __m128 squared_dot_product(const float* v1, const float* v2, OutputControl control) noexcept {
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
                }
            #else
                sv1[3] = sv2[3] = 0;                        // zero out the weights
                __m128 diff = _mm_sub_ps(sv1, sv2);         // calculate the difference
                __m128 multiplied = _mm_mul_ps(diff, diff); // square the difference
                return _mm_hadd_ps(multiplied, multiplied); // sum the components
            #endif
        }
    }

    inline ausaxs::hist::detail::EvaluatedResult ausaxs::hist::detail::CompactCoordinatesData::evaluate_sse(
        const CompactCoordinatesData& other
    ) const noexcept {
        __m128 dist2 = squared_dot_product(this->data.data(), other.data.data(), OutputControl::ALL);
        __m128 dist_sqrt = _mm_sqrt_ps(dist2);
        float dist = _mm_cvtss_f32(dist_sqrt);
        return EvaluatedResult(dist, this->value.w*other.value.w);
    }

    inline ausaxs::hist::detail::EvaluatedResultRounded ausaxs::hist::detail::CompactCoordinatesData::evaluate_rounded_sse(
        const CompactCoordinatesData& other
    ) const noexcept {
        __m128 dist2 = squared_dot_product(this->data.data(), other.data.data(), OutputControl::ALL);
        __m128 dist_sqrt = _mm_sqrt_ps(dist2);
        int32_t dist_bin = std::round(get_inv_width()*_mm_cvtss_f32(dist_sqrt));
        return EvaluatedResultRounded(dist_bin, this->value.w*other.value.w);
    }

    template<bool vbw>
    inline ausaxs::hist::detail::QuadEvaluatedResult ausaxs::hist::detail::CompactCoordinatesData<vbw>::evaluate_sse(
        const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4
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

        QuadEvaluatedResult result;
        _mm_store_ps(reinterpret_cast<float*>(result.distances.data()), dist_sqrt);         // efficient store of distances
        _mm_store_ps(reinterpret_cast<float*>(result.weights.data()), weights);             // efficient store of weights
        return result;
    }

    template<bool vbw>
    inline ausaxs::hist::detail::QuadEvaluatedResultRounded ausaxs::hist::detail::CompactCoordinatesData<vbw>::evaluate_rounded_sse(
        const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4
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

        QuadEvaluatedResultRounded result;
        _mm_store_si128(reinterpret_cast<__m128i*>(result.distances.data()), dist_bin); // efficient store of bins
        _mm_store_ps(reinterpret_cast<float*>(result.weights.data()), weights);         // efficient store of weights
        return result;
    }

    template<bool vbw>
    inline ausaxs::hist::detail::OctoEvaluatedResult ausaxs::hist::detail::CompactCoordinatesData<vbw>::evaluate_sse(
        const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4, 
        const CompactCoordinatesData& v5, const CompactCoordinatesData& v6, const CompactCoordinatesData& v7, const CompactCoordinatesData& v8
    ) const noexcept {
        OctoEvaluatedResult result;
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
    inline ausaxs::hist::detail::OctoEvaluatedResultRounded ausaxs::hist::detail::CompactCoordinatesData<vbw>::evaluate_rounded_sse(
        const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4, 
        const CompactCoordinatesData& v5, const CompactCoordinatesData& v6, const CompactCoordinatesData& v7, const CompactCoordinatesData& v8
    ) const noexcept {
        OctoEvaluatedResultRounded result;
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
    #include <immintrin.h>

    namespace ausaxs::hist::detail {
        /**
        * @brief Calculate the squared distance between three CompactCoordinatesData using AVX instructions.
        */
        static inline __m256 squared_dot_product(const float* v, const float* v1, const float* v2, OutputControl control) noexcept {
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
            }
        }
    }

    inline ausaxs::hist::detail::EvaluatedResult ausaxs::hist::detail::CompactCoordinatesData::evaluate_avx(
        const CompactCoordinatesData& other
    ) const noexcept {
        return evaluate_sse(other); // no way to optimize a single evaluation with AVX
    }

    inline ausaxs::hist::detail::EvaluatedResultRounded ausaxs::hist::detail::CompactCoordinatesData::evaluate_rounded_avx(
        const CompactCoordinatesData& other
    ) const noexcept {
        return evaluate_rounded_sse(other); // no way to optimize a single evaluation with AVX
    }

    template<bool vbw>
    inline ausaxs::hist::detail::QuadEvaluatedResult ausaxs::hist::detail::CompactCoordinatesData<vbw>::evaluate_avx(
        const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4
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

        QuadEvaluatedResult result;
        _mm_store_ps(reinterpret_cast<float*>(result.distances.data()), dist_sqrt);         // efficient store of distances
        _mm_store_ps(reinterpret_cast<float*>(result.weights.data()), weights);             // efficient store of weights
        return result;
    }

    template<bool vbw>
    inline ausaxs::hist::detail::QuadEvaluatedResultRounded ausaxs::hist::detail::CompactCoordinatesData<vbw>::evaluate_rounded_avx(
        const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4
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

        QuadEvaluatedResultRounded result;
        _mm_store_si128(reinterpret_cast<__m128i*>(result.distances.data()), dist_bin);  // efficient store of bins
        _mm_store_ps(reinterpret_cast<float*>(result.weights.data()), weights);          // efficient store of weights
        return result;
    }

    template<bool vbw>
    inline ausaxs::hist::detail::OctoEvaluatedResult ausaxs::hist::detail::CompactCoordinatesData<vbw>::evaluate_avx(
        const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4, 
        const CompactCoordinatesData& v5, const CompactCoordinatesData& v6, const CompactCoordinatesData& v7, const CompactCoordinatesData& v8
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

        OctoEvaluatedResult result;
        _mm256_store_ps(reinterpret_cast<float*>(result.distances.data()), dist_sqrt);          // efficient store of distances
        _mm256_store_ps(reinterpret_cast<float*>(result.weights.data()), weights);              // efficient store of weights
        return result;
    }

    template<bool vbw>
    inline ausaxs::hist::detail::OctoEvaluatedResultRounded ausaxs::hist::detail::CompactCoordinatesData<vbw>::evaluate_rounded_avx(
        const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4, 
        const CompactCoordinatesData& v5, const CompactCoordinatesData& v6, const CompactCoordinatesData& v7, const CompactCoordinatesData& v8
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

        OctoEvaluatedResultRounded result;
        _mm256_store_si256(reinterpret_cast<__m256i*>(result.distances.data()), dist_bin);   // efficient store of bins
        _mm256_store_ps(reinterpret_cast<float*>(result.weights.data()), weights);           // efficient store of weights
        return result;
    }
#endif