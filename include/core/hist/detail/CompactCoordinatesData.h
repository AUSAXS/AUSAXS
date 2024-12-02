#pragma once

#include <math/Vector3.h>

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
        EvaluatedResult() = default;
        EvaluatedResult(float distance, float weight) : distance(distance), weight(weight) {}
        float distance;         // The raw distance 
        float weight;           // The combined weight
    };

    struct EvaluatedResultRounded {
        EvaluatedResultRounded() = default;
        EvaluatedResultRounded(int32_t distance_bin, float weight) : distance(distance_bin), weight(weight) {}
        int32_t distance;       // The distance bin 
        float weight;           // The combined weight
    };

    /**
     * @brief Simple structure for storing the results of four distance and weight calculations.
     *        This is necessary to restructure the memory layout for more efficient 128-bit SIMD instructions.
     */
    struct QuadEvaluatedResult {
        QuadEvaluatedResult() = default;
        QuadEvaluatedResult(const EvaluatedResult& v1, const EvaluatedResult& v2, const EvaluatedResult& v3, const EvaluatedResult& v4) 
            : distances{v1.distance, v2.distance, v3.distance, v4.distance}, weights{v1.weight, v2.weight, v3.weight, v4.weight}
        {}
        QuadEvaluatedResult(const std::array<float, 4>& distances, const std::array<float, 4>& weights) 
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
        QuadEvaluatedResultRounded() = default;
        QuadEvaluatedResultRounded(const EvaluatedResultRounded& v1, const EvaluatedResultRounded& v2, const EvaluatedResultRounded& v3, const EvaluatedResultRounded& v4) 
            : distances{v1.distance, v2.distance, v3.distance, v4.distance}, 
                weights{v1.weight, v2.weight, v3.weight, v4.weight} 
        {}
        QuadEvaluatedResultRounded(const std::array<int32_t, 4>& distances, const std::array<float, 4>& weights) 
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
        OctoEvaluatedResult() = default;
        OctoEvaluatedResult(
            const EvaluatedResult& v1, const EvaluatedResult& v2, const EvaluatedResult& v3, const EvaluatedResult& v4, 
            const EvaluatedResult& v5, const EvaluatedResult& v6, const EvaluatedResult& v7, const EvaluatedResult& v8) 
            : distances{v1.distance, v2.distance, v3.distance, v4.distance, v5.distance, v6.distance, v7.distance, v8.distance},
              weights{v1.weight, v2.weight, v3.weight, v4.weight, v5.weight, v6.weight, v7.weight, v8.weight}
        {}
        OctoEvaluatedResult(const std::array<float, 8>& distances, const std::array<float, 8>& weights) 
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
        OctoEvaluatedResultRounded() = default;
        OctoEvaluatedResultRounded(
            const EvaluatedResultRounded& v1, const EvaluatedResultRounded& v2, const EvaluatedResultRounded& v3, const EvaluatedResultRounded& v4, 
            const EvaluatedResultRounded& v5, const EvaluatedResultRounded& v6, const EvaluatedResultRounded& v7, const EvaluatedResultRounded& v8) 
            : distances{v1.distance, v2.distance, v3.distance, v4.distance, v5.distance, v6.distance, v7.distance, v8.distance}, 
              weights{v1.weight, v2.weight, v3.weight, v4.weight, v5.weight, v6.weight, v7.weight, v8.weight} 
        {}
        OctoEvaluatedResultRounded(const std::array<int32_t, 8>& distances, const std::array<float, 8>& weights) 
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

    class CompactCoordinatesData {
        public:
            CompactCoordinatesData() = default;

            template<numeric T, numeric V>
            CompactCoordinatesData(const Vector3<T>& v, V w) : value{.pos={static_cast<float>(v.x()), static_cast<float>(v.y()), static_cast<float>(v.z())}, .w=static_cast<float>(w)} {}
            CompactCoordinatesData(const Vector3<float>& v, float w) : value{.pos=v, .w=w} {}

            /**
             * @brief Calculate the @a binned distance and combined weight between this and a single other CompactCoordinatesData.
             */
            EvaluatedResultRounded evaluate_rounded(const CompactCoordinatesData& other) const;

            /**
             * @brief Calculate the distance and combined weight between this and a single other CompactCoordinatesData.
             */
            EvaluatedResult evaluate(const CompactCoordinatesData& other) const;

            /**
             * @brief Calculate the @a binned distance and combined weight between this and four other CompactCoordinatesData.
             *        This leverages more efficient SIMD instructions by using a 128-bit registers (SSE).
             */
            QuadEvaluatedResultRounded evaluate_rounded(const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4) const;

            /**
             * @brief Calculate the distance and combined weight between this and four other CompactCoordinatesData.
             *        This leverages more efficient SIMD instructions by using a 128-bit registers (SSE).
             */
            QuadEvaluatedResult evaluate(const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4) const;

            /**
             * @brief Calculate the @a binned distance and combined weight between this and four other CompactCoordinatesData.
             *        This leverages more efficient SIMD instructions by using either two 128-bit registers (SSE) or one 256-bit register (AVX).
             */
            OctoEvaluatedResultRounded evaluate_rounded(
                const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4,
                const CompactCoordinatesData& v5, const CompactCoordinatesData& v6, const CompactCoordinatesData& v7, const CompactCoordinatesData& v8
            ) const;

            /**
             * @brief Calculate the distance and combined weight between this and four other CompactCoordinatesData.
             *        This leverages more efficient SIMD instructions by using either two 128-bit registers (SSE) or one 256-bit register (AVX).
             */
            OctoEvaluatedResult evaluate(
                const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4,
                const CompactCoordinatesData& v5, const CompactCoordinatesData& v6, const CompactCoordinatesData& v7, const CompactCoordinatesData& v8
            ) const;

            union {
                struct {Vector3<float> pos; float w;} value;
                std::array<float, 4> data;
            };

        protected:
            EvaluatedResultRounded evaluate_rounded_scalar(const CompactCoordinatesData& other) const;
            EvaluatedResult evaluate_scalar(const CompactCoordinatesData& other) const;

            QuadEvaluatedResultRounded evaluate_rounded_scalar(const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4) const;
            QuadEvaluatedResult evaluate_scalar(const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4) const;

            OctoEvaluatedResultRounded evaluate_rounded_scalar(
                const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4,
                const CompactCoordinatesData& v5, const CompactCoordinatesData& v6, const CompactCoordinatesData& v7, const CompactCoordinatesData& v8
            ) const;
            OctoEvaluatedResult evaluate_scalar(
                const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4,
                const CompactCoordinatesData& v5, const CompactCoordinatesData& v6, const CompactCoordinatesData& v7, const CompactCoordinatesData& v8
            ) const;

            #if defined __SSE2__
                EvaluatedResultRounded evaluate_rounded_sse(const CompactCoordinatesData& other) const;
                EvaluatedResult evaluate_sse(const CompactCoordinatesData& other) const;

                QuadEvaluatedResultRounded evaluate_rounded_sse(const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4) const;
                QuadEvaluatedResult evaluate_sse(const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4) const;

                OctoEvaluatedResultRounded evaluate_rounded_sse(
                    const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4,
                    const CompactCoordinatesData& v5, const CompactCoordinatesData& v6, const CompactCoordinatesData& v7, const CompactCoordinatesData& v8
                ) const;
                OctoEvaluatedResult evaluate_sse(
                    const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4,
                    const CompactCoordinatesData& v5, const CompactCoordinatesData& v6, const CompactCoordinatesData& v7, const CompactCoordinatesData& v8
                ) const;
            #endif

            #if defined __AVX__
                EvaluatedResultRounded evaluate_rounded_avx(const CompactCoordinatesData& other) const;
                EvaluatedResult evaluate_avx(const CompactCoordinatesData& other) const;

                QuadEvaluatedResultRounded evaluate_rounded_avx(const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4) const;
                QuadEvaluatedResult evaluate_avx(const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4) const;

                OctoEvaluatedResultRounded evaluate_rounded_avx(
                    const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4,
                    const CompactCoordinatesData& v5, const CompactCoordinatesData& v6, const CompactCoordinatesData& v7, const CompactCoordinatesData& v8
                ) const;
                OctoEvaluatedResult evaluate_avx(
                    const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4,
                    const CompactCoordinatesData& v5, const CompactCoordinatesData& v6, const CompactCoordinatesData& v7, const CompactCoordinatesData& v8
                ) const;
            #endif
    };
    static_assert(sizeof(CompactCoordinatesData) == 16,              "CompactCoordinatesData is not 16 bytes. This is required for aligning SIMD instructions.");
    static_assert(std::is_trivial_v<CompactCoordinatesData>,         "CompactCoordinatesData is not trivial");
    static_assert(std::is_standard_layout_v<CompactCoordinatesData>, "CompactCoordinatesData is not standard layout");
    static_assert(supports_nothrow_move_v<CompactCoordinatesData>,   "CompactCoordinatesData should support nothrow move semantics.");
}