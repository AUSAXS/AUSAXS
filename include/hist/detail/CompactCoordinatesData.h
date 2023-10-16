#pragma once

#include <math/MathFwd.h>

#include <array>
#include <cstdint>

namespace hist::detail {
    /**
     * @brief Simple structure for storing the results of a distance and weight calculation.
     */
    struct EvaluatedResult {
        EvaluatedResult() = default;
        EvaluatedResult(int32_t distance, float weight) : distance(distance), weight(weight) {}
        int32_t distance; // The distance bin 
        float weight;     // The combined weight
    };

    /**
     * @brief Simple structure for storing the results of four distance and weight calculations.
     *        This is necessary to restructure the memory layout for more efficient 128-bit SIMD instructions.
     */
    struct QuadEvaluatedResult {
        QuadEvaluatedResult() = default;
        QuadEvaluatedResult(const EvaluatedResult& v1, const EvaluatedResult& v2, const EvaluatedResult& v3, const EvaluatedResult& v4) : distance{v1.distance, v2.distance, v3.distance, v4.distance}, weight{v1.weight, v2.weight, v3.weight, v4.weight} {}
        union {
            struct {int32_t first, second, third, fourth;} distances; // The distance bin
            std::array<int32_t, 4> distance;                          // The distance bin
        };
        union {
            struct {float first, second, third, fourth;} weights;   // The combined weight
            std::array<float, 4> weight;                            // The combined weight
        };
    };

    /**
     * @brief Simple structure for storing the results of eight distance and weight calculations.
     *        This is necessary to restructure the memory layout for more efficient 256-bit SIMD instructions.
     */
    struct OctoEvaluatedResult {
        OctoEvaluatedResult() = default;
        OctoEvaluatedResult(const EvaluatedResult& v1, const EvaluatedResult& v2, const EvaluatedResult& v3, const EvaluatedResult& v4, const EvaluatedResult& v5, const EvaluatedResult& v6, const EvaluatedResult& v7, const EvaluatedResult& v8) : distance{v1.distance, v2.distance, v3.distance, v4.distance, v5.distance, v6.distance, v7.distance, v8.distance}, weight{v1.weight, v2.weight, v3.weight, v4.weight, v5.weight, v6.weight, v7.weight, v8.weight} {}
        union {
            struct {int32_t first, second, third, fourth, fifth, sixth, seventh, eighth;} distances; // The distance bin
            std::array<int32_t, 8> distance;                                                         // The distance bin
        };
        union {
            struct {float first, second, third, fourth, fifth, sixth, seventh, eighth;} weights; // The combined weight
            std::array<float, 8> weight;                                                         // The combined weight
        };
    };

    class CompactCoordinatesData {
        public:
            CompactCoordinatesData();

            CompactCoordinatesData(const Vector3<double>& v, float w);

            /**
             * @brief Calculate the distance and combined weight between this and a single other CompactCoordinatesData.
             */
            EvaluatedResult evaluate(const CompactCoordinatesData& other) const;

            /**
             * @brief Calculate the distance and combined weight between this and four other CompactCoordinatesData.
             *        This leverages more efficient SIMD instructions by using a 128-bit registers (SSE).
             */
            QuadEvaluatedResult evaluate(const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4) const;

            /**
             * @brief Calculate the distance and combined weight between this and four other CompactCoordinatesData.
             *        This leverages more efficient SIMD instructions by using either two 128-bit registers (SSE) or one 256-bit register (AVX).
             */
            OctoEvaluatedResult evaluate(const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4,
                                        const CompactCoordinatesData& v5, const CompactCoordinatesData& v6, const CompactCoordinatesData& v7, const CompactCoordinatesData& v8) const;

            union {
                struct {float x, y, z, w;} value;
                std::array<float, 4> data;
            };

        protected:
            EvaluatedResult evaluate_scalar(const CompactCoordinatesData& other) const;
            QuadEvaluatedResult evaluate_scalar(const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4) const;
            OctoEvaluatedResult evaluate_scalar(const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4,
                                                const CompactCoordinatesData& v5, const CompactCoordinatesData& v6, const CompactCoordinatesData& v7, const CompactCoordinatesData& v8) const;

            #if defined __SSE2__
                EvaluatedResult evaluate_sse(const CompactCoordinatesData& other) const;
                QuadEvaluatedResult evaluate_sse(const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4) const;
                OctoEvaluatedResult evaluate_sse(const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4,
                                                const CompactCoordinatesData& v5, const CompactCoordinatesData& v6, const CompactCoordinatesData& v7, const CompactCoordinatesData& v8) const;
            #endif

            #if defined __AVX__
                EvaluatedResult evaluate_avx(const CompactCoordinatesData& other) const;
                QuadEvaluatedResult evaluate_avx(const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4) const;
                OctoEvaluatedResult evaluate_avx(const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4,
                                                const CompactCoordinatesData& v5, const CompactCoordinatesData& v6, const CompactCoordinatesData& v7, const CompactCoordinatesData& v8) const;            
            #endif
    };
    static_assert(sizeof(CompactCoordinatesData) == 16, "hist::detail::CompactCoordinates::CompactCoordinatesData is not 16 bytes. This is required for aligning SIMD instructions.");
}