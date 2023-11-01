#include <hist/detail/CompactCoordinatesData.h>
#include <constants/Constants.h>
#include <math/Vector3.h>

using namespace hist::detail;

constexpr float inv_width = 1.0f/constants::axes::d_axis.width();

CompactCoordinatesData::CompactCoordinatesData() : data(std::array<float, 4>({0, 0, 0, 0})) {}

EvaluatedResult CompactCoordinatesData::evaluate(const CompactCoordinatesData& other) const {
    #if defined __SSE2__
        return evaluate_sse(other);
    #else
        return evaluate_scalar(other);
    #endif
}

EvaluatedResultRounded CompactCoordinatesData::evaluate_rounded(const CompactCoordinatesData& other) const {
    #if defined __SSE2__
        return evaluate_rounded_sse(other);
    #else
        return evaluate_rounded_scalar(other);
    #endif
}

QuadEvaluatedResult CompactCoordinatesData::evaluate(const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4) const {
    #if defined __AVX__
        return evaluate_avx(v1, v2, v3, v4);
    #elif defined __SSE2__
        return evaluate_sse(v1, v2, v3, v4);
    #else
        return evaluate_scalar(v1, v2, v3, v4);
    #endif
}

QuadEvaluatedResultRounded CompactCoordinatesData::evaluate_rounded(const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4) const {
    #if defined __AVX__
        return evaluate_rounded_avx(v1, v2, v3, v4);
    #elif defined __SSE2__
        return evaluate_rounded_sse(v1, v2, v3, v4);
    #else
        return evaluate_rounded_scalar(v1, v2, v3, v4);
    #endif
}

OctoEvaluatedResult CompactCoordinatesData::evaluate(const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4, const CompactCoordinatesData& v5, const CompactCoordinatesData& v6, const CompactCoordinatesData& v7, const CompactCoordinatesData& v8) const {
    #if defined __AVX__
        return evaluate_avx(v1, v2, v3, v4, v5, v6, v7, v8);
    #elif defined __SSE2__
        return evaluate_sse(v1, v2, v3, v4, v5, v6, v7, v8);
    #else
        return evaluate_scalar(v1, v2, v3, v4, v5, v6, v7, v8);
    #endif
}

OctoEvaluatedResultRounded CompactCoordinatesData::evaluate_rounded(const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4, const CompactCoordinatesData& v5, const CompactCoordinatesData& v6, const CompactCoordinatesData& v7, const CompactCoordinatesData& v8) const {
    #if defined __AVX__
        return evaluate_rounded_avx(v1, v2, v3, v4, v5, v6, v7, v8);
    #elif defined __SSE2__
        return evaluate_rounded_sse(v1, v2, v3, v4, v5, v6, v7, v8);
    #else
        return evaluate_rounded_scalar(v1, v2, v3, v4, v5, v6, v7, v8);
    #endif
}

static inline float squared_dot_product(const float* v1, const float* v2) {
    float dx = v1[0] - v2[0];
    float dy = v1[1] - v2[1];
    float dz = v1[2] - v2[2];
    return dx*dx + dy*dy + dz*dz;
}

EvaluatedResult CompactCoordinatesData::evaluate_scalar(const CompactCoordinatesData& other) const {
    float dist = std::sqrt(squared_dot_product(this->data.data(), other.data.data()));
    return EvaluatedResult(dist, value.w*other.value.w);
}

EvaluatedResultRounded CompactCoordinatesData::evaluate_rounded_scalar(const CompactCoordinatesData& other) const {
    int32_t dist = std::round(inv_width*std::sqrt(squared_dot_product(this->data.data(), other.data.data())));
    return EvaluatedResultRounded(dist, value.w*other.value.w);
}

QuadEvaluatedResult CompactCoordinatesData::evaluate_scalar(const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4) const {
    float dx1 = std::sqrt(squared_dot_product(this->data.data(), v1.data.data()));
    float dx2 = std::sqrt(squared_dot_product(this->data.data(), v2.data.data()));
    float dx3 = std::sqrt(squared_dot_product(this->data.data(), v3.data.data()));
    float dx4 = std::sqrt(squared_dot_product(this->data.data(), v4.data.data()));
    return QuadEvaluatedResult({dx1, value.w*v1.value.w}, {dx2, value.w*v2.value.w}, {dx3, value.w*v3.value.w}, {dx4, value.w*v4.value.w});
}

QuadEvaluatedResultRounded CompactCoordinatesData::evaluate_rounded_scalar(const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4) const {
    int32_t dx1 = std::round(inv_width*std::sqrt(squared_dot_product(this->data.data(), v1.data.data())));
    int32_t dx2 = std::round(inv_width*std::sqrt(squared_dot_product(this->data.data(), v2.data.data())));
    int32_t dx3 = std::round(inv_width*std::sqrt(squared_dot_product(this->data.data(), v3.data.data())));
    int32_t dx4 = std::round(inv_width*std::sqrt(squared_dot_product(this->data.data(), v4.data.data())));
    return QuadEvaluatedResultRounded({dx1, value.w*v1.value.w}, {dx2, value.w*v2.value.w}, {dx3, value.w*v3.value.w}, {dx4, value.w*v4.value.w});
}

OctoEvaluatedResult CompactCoordinatesData::evaluate_scalar(const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4, const CompactCoordinatesData& v5, const CompactCoordinatesData& v6, const CompactCoordinatesData& v7, const CompactCoordinatesData& v8) const {
    float dx1 = std::sqrt(squared_dot_product(this->data.data(), v1.data.data()));
    float dx2 = std::sqrt(squared_dot_product(this->data.data(), v2.data.data()));
    float dx3 = std::sqrt(squared_dot_product(this->data.data(), v3.data.data()));
    float dx4 = std::sqrt(squared_dot_product(this->data.data(), v4.data.data()));
    float dx5 = std::sqrt(squared_dot_product(this->data.data(), v5.data.data()));
    float dx6 = std::sqrt(squared_dot_product(this->data.data(), v6.data.data()));
    float dx7 = std::sqrt(squared_dot_product(this->data.data(), v7.data.data()));
    float dx8 = std::sqrt(squared_dot_product(this->data.data(), v8.data.data()));
    return OctoEvaluatedResult({dx1, value.w*v1.value.w}, {dx2, value.w*v2.value.w}, {dx3, value.w*v3.value.w}, {dx4, value.w*v4.value.w}, {dx5, value.w*v5.value.w}, {dx6, value.w*v6.value.w}, {dx7, value.w*v7.value.w}, {dx8, value.w*v8.value.w});
}

OctoEvaluatedResultRounded CompactCoordinatesData::evaluate_rounded_scalar(const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4, const CompactCoordinatesData& v5, const CompactCoordinatesData& v6, const CompactCoordinatesData& v7, const CompactCoordinatesData& v8) const {
    int32_t dx1 = std::round(inv_width*std::sqrt(squared_dot_product(this->data.data(), v1.data.data())));
    int32_t dx2 = std::round(inv_width*std::sqrt(squared_dot_product(this->data.data(), v2.data.data())));
    int32_t dx3 = std::round(inv_width*std::sqrt(squared_dot_product(this->data.data(), v3.data.data())));
    int32_t dx4 = std::round(inv_width*std::sqrt(squared_dot_product(this->data.data(), v4.data.data())));
    int32_t dx5 = std::round(inv_width*std::sqrt(squared_dot_product(this->data.data(), v5.data.data())));
    int32_t dx6 = std::round(inv_width*std::sqrt(squared_dot_product(this->data.data(), v6.data.data())));
    int32_t dx7 = std::round(inv_width*std::sqrt(squared_dot_product(this->data.data(), v7.data.data())));
    int32_t dx8 = std::round(inv_width*std::sqrt(squared_dot_product(this->data.data(), v8.data.data())));
    return OctoEvaluatedResultRounded({dx1, value.w*v1.value.w}, {dx2, value.w*v2.value.w}, {dx3, value.w*v3.value.w}, {dx4, value.w*v4.value.w}, {dx5, value.w*v5.value.w}, {dx6, value.w*v6.value.w}, {dx7, value.w*v7.value.w}, {dx8, value.w*v8.value.w});
}

#if defined __SSE2__
    #include <nmmintrin.h>
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
    static inline __m128 squared_dot_product(const float* v1, const float* v2, OutputControl control) {
        // load data into SSE registers
        __m128 sv1 = _mm_load_ps(v1);
        __m128 sv2 = _mm_load_ps(v2);
        #if defined __SSE4_1__
            __m128 diff = _mm_sub_ps(sv1, sv2);         // calculate the difference
            return _mm_dp_ps(diff, diff, control);      // calculate the dot product, masked to only use the first three components
        #else
            sv1[3] = sv2[3] = 0;                        // zero out the weights
            __m128 diff = _mm_sub_ps(sv1, sv2);         // calculate the difference
            __m128 multiplied = _mm_mul_ps(diff, diff); // square the difference
            return _mm_hadd_ps(multiplied, multiplied); // sum the components
        #endif
    }

    EvaluatedResult CompactCoordinatesData::evaluate_sse(const CompactCoordinatesData& other) const {
        __m128 dist2 = squared_dot_product(this->data.data(), other.data.data(), OutputControl::ALL);
        __m128 dist_sqrt = _mm_sqrt_ps(dist2);
        float dist = _mm_cvtss_f32(dist_sqrt);
        return EvaluatedResult(dist, this->value.w*other.value.w);
    }

    EvaluatedResultRounded CompactCoordinatesData::evaluate_rounded_sse(const CompactCoordinatesData& other) const {
        auto result = evaluate_sse(other);
        int32_t dist_bin = std::round(inv_width*result.distance);
        return EvaluatedResultRounded(dist_bin, result.weight);
    }

    QuadEvaluatedResult CompactCoordinatesData::evaluate_sse(const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4) const {
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
        _mm_store_ps(reinterpret_cast<float*>(result.distance.data()), dist_sqrt); // efficient store of distances
        _mm_store_ps(reinterpret_cast<float*>(result.weight.data()), weights);     // efficient store of weights
        return result;
    }

    QuadEvaluatedResultRounded CompactCoordinatesData::evaluate_rounded_sse(const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4) const {
        auto result = evaluate_sse(v1, v2, v3, v4);
        __m128 dist_binf = _mm_mul_ps(*reinterpret_cast<__m128*>(result.distance.data()), _mm_set_ps1(inv_width));  // multiply by the inverse width
        __m128i dist_bin = _mm_cvtps_epi32(dist_binf);                                                              // cast to int

        _mm_store_si128(reinterpret_cast<__m128i*>(result.distance.data()), dist_bin); // efficient store of bins

        return reinterpret_cast<QuadEvaluatedResultRounded&>(result);
    }

    OctoEvaluatedResult CompactCoordinatesData::evaluate_sse(const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4, const CompactCoordinatesData& v5, const CompactCoordinatesData& v6, const CompactCoordinatesData& v7, const CompactCoordinatesData& v8) const {
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

            _mm_store_ps(reinterpret_cast<float*>(result.distance.data()), dist_sqrt);
            _mm_store_ps(reinterpret_cast<float*>(result.weight.data()), weights);
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

            _mm_store_ps(reinterpret_cast<float*>(result.distance.data()+4), dist_sqrt);
            _mm_store_ps(reinterpret_cast<float*>(result.weight.data()+4), weights);
        }
        return result;
    }

    OctoEvaluatedResultRounded CompactCoordinatesData::evaluate_rounded_sse(const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4, const CompactCoordinatesData& v5, const CompactCoordinatesData& v6, const CompactCoordinatesData& v7, const CompactCoordinatesData& v8) const {
        auto result = evaluate_sse(v1, v2, v3, v4, v5, v6, v7, v8);
        {   // first four
            __m128 dist_binf = _mm_mul_ps(*reinterpret_cast<__m128*>(result.distance.data()), _mm_set_ps1(inv_width));  // multiply by the inverse width
            __m128i dist_bin = _mm_cvtps_epi32(dist_binf);                                                              // cast to int
            _mm_store_si128(reinterpret_cast<__m128i*>(result.distance.data()), dist_bin);
        }
        {   // last four
            __m128 dist_binf = _mm_mul_ps(*reinterpret_cast<__m128*>(result.distance.data()+4), _mm_set_ps1(inv_width));  // multiply by the inverse width
            __m128i dist_bin = _mm_cvtps_epi32(dist_binf);                                                                // cast to int
            _mm_store_si128(reinterpret_cast<__m128i*>(result.distance.data()+4), dist_bin);
        }
        return reinterpret_cast<OctoEvaluatedResultRounded&>(result);
    }
#endif

#if defined __AVX__
    #include <immintrin.h>

    /**
     * @brief Calculate the squared distance between three CompactCoordinatesData using AVX instructions.
     */
    static inline __m256 squared_dot_product(const float* v, const float* v1, const float* v2, OutputControl control) {
        // load data into the 256 bit registers
        __m128 sv = _mm_load_ps(v);
        __m128 sv1 = _mm_load_ps(v1);
        __m128 sv2 = _mm_load_ps(v2);

        __m256 svv = _mm256_broadcast_ps(&sv);      // copy the first 128 bits to the second 128 bits
        __m256 sv12 = _mm256_set_m128(sv2, sv1);    // combine the two 128 bit registers
        __m256 diff = _mm256_sub_ps(svv, sv12);     // calculate the difference
        return _mm256_dp_ps(diff, diff, control);   // calculate the dot product, masked to only use the first three components
    }

    EvaluatedResult CompactCoordinatesData::evaluate_avx(const CompactCoordinatesData& other) const {
        return evaluate_sse(other); // no way to optimize a single evaluation with AVX
    }

    EvaluatedResultRounded CompactCoordinatesData::evaluate_rounded_avx(const CompactCoordinatesData& other) const {
        return evaluate_rounded_sse(other); // no way to optimize a single evaluation with AVX
    }

    QuadEvaluatedResult CompactCoordinatesData::evaluate_avx(const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4) const {
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
        _mm_store_ps(reinterpret_cast<float*>(result.distance.data()), dist_sqrt); // efficient store of bins
        _mm_store_ps(reinterpret_cast<float*>(result.weight.data()), weights);     // efficient store of weights
        return result;
    }

    QuadEvaluatedResultRounded CompactCoordinatesData::evaluate_rounded_avx(const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4) const {
        auto result = evaluate_avx(v1, v2, v3, v4);
        __m128 dist_binf = _mm_mul_ps(*reinterpret_cast<__m128*>(result.distance.data()), _mm_set_ps1(inv_width));  // multiply by the inverse width
        __m128i dist_bin = _mm_cvtps_epi32(dist_binf);                                                              // cast to int

        _mm_store_si128(reinterpret_cast<__m128i*>(result.distance.data()), dist_bin); // efficient store of bins
        return reinterpret_cast<QuadEvaluatedResultRounded&>(result);
    }

    OctoEvaluatedResult CompactCoordinatesData::evaluate_avx(const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4, const CompactCoordinatesData& v5, const CompactCoordinatesData& v6, const CompactCoordinatesData& v7, const CompactCoordinatesData& v8) const {
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
        _mm256_store_ps(reinterpret_cast<float*>(result.distance.data()), dist_sqrt); // efficient store of bins
        _mm256_store_ps(reinterpret_cast<float*>(result.weight.data()), weights);     // efficient store of weights
        return result;
    }

    OctoEvaluatedResultRounded CompactCoordinatesData::evaluate_rounded_avx(const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4, const CompactCoordinatesData& v5, const CompactCoordinatesData& v6, const CompactCoordinatesData& v7, const CompactCoordinatesData& v8) const {
        auto result = evaluate_avx(v1, v2, v3, v4, v5, v6, v7, v8);
        __m256 dist_binf = _mm256_mul_ps(*reinterpret_cast<__m256*>(result.distance.data()), _mm256_set1_ps(inv_width)); // multiply by the inverse width
        __m256i dist_bin = _mm256_cvtps_epi32(dist_binf);                                                                // cast to int

        _mm256_store_si256(reinterpret_cast<__m256i*>(result.distance.data()), dist_bin); // efficient store of bins
        return reinterpret_cast<OctoEvaluatedResultRounded&>(result);
    }
#endif