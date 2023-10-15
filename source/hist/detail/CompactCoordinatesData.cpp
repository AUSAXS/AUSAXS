#include <hist/detail/CompactCoordinatesData.h>
#include <math/Vector3.h>

using namespace hist::detail;

CompactCoordinatesData::CompactCoordinatesData() : data(std::array<float, 4>(0)) {}

CompactCoordinatesData::CompactCoordinatesData(const Vector3<double>& v, float w) : value{.x=v.x(), .y=v.y(), .z=v.z(), .w=w} {}

EvaluatedResult CompactCoordinatesData::evaluate(const CompactCoordinatesData& other) const {
    #if defined __SSE2__
        return evaluate_sse(other);
    #else
        return evaluate_scalar(other);
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

OctoEvaluatedResult CompactCoordinatesData::evaluate(const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4, const CompactCoordinatesData& v5, const CompactCoordinatesData& v6, const CompactCoordinatesData& v7, const CompactCoordinatesData& v8) const {
    #if defined __AVX__
        return evaluate_avx(v1, v2, v3, v4, v5, v6, v7, v8);
    #elif defined __SSE2__
        return evaluate_sse(v1, v2, v3, v4, v5, v6, v7, v8);
    #else
        return evaluate_scalar(v1, v2, v3, v4, v5, v6, v7, v8);
    #endif
}

static inline float squared_dot_product(const float* v1, const float* v2) {
    float dx = v1[0] - v2[0];
    float dy = v1[1] - v2[1];
    float dz = v1[2] - v2[2];
    return dx*dx + dy*dy + dz*dz;
}

EvaluatedResult CompactCoordinatesData::evaluate_scalar(const CompactCoordinatesData& other) const {
    float dist = inv_width*std::sqrt(squared_dot_product(this->data.data(), other.data.data()));
    return EvaluatedResult(dist, value.w*other.value.w);
}

QuadEvaluatedResult CompactCoordinatesData::evaluate_scalar(const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4) const {
    float dx1 = inv_width*std::sqrt(squared_dot_product(this->data.data(), v1.data.data()));
    float dx2 = inv_width*std::sqrt(squared_dot_product(this->data.data(), v2.data.data()));
    float dx3 = inv_width*std::sqrt(squared_dot_product(this->data.data(), v3.data.data()));
    float dx4 = inv_width*std::sqrt(squared_dot_product(this->data.data(), v4.data.data()));
    return QuadEvaluatedResult({dx1, value.w*v1.value.w}, {dx2, value.w*v2.value.w}, {dx3, value.w*v3.value.w}, {dx4, value.w*v4.value.w});
}
OctoEvaluatedResult CompactCoordinatesData::evaluate_scalar(const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4, const CompactCoordinatesData& v5, const CompactCoordinatesData& v6, const CompactCoordinatesData& v7, const CompactCoordinatesData& v8) const {
    float dx1 = inv_width*std::sqrt(squared_dot_product(this->data.data(), v1.data.data()));
    float dx2 = inv_width*std::sqrt(squared_dot_product(this->data.data(), v2.data.data()));
    float dx3 = inv_width*std::sqrt(squared_dot_product(this->data.data(), v3.data.data()));
    float dx4 = inv_width*std::sqrt(squared_dot_product(this->data.data(), v4.data.data()));
    float dx5 = inv_width*std::sqrt(squared_dot_product(this->data.data(), v5.data.data()));
    float dx6 = inv_width*std::sqrt(squared_dot_product(this->data.data(), v6.data.data()));
    float dx7 = inv_width*std::sqrt(squared_dot_product(this->data.data(), v7.data.data()));
    float dx8 = inv_width*std::sqrt(squared_dot_product(this->data.data(), v8.data.data()));
    return OctoEvaluatedResult({dx1, value.w*v1.value.w}, {dx2, value.w*v2.value.w}, {dx3, value.w*v3.value.w}, {dx4, value.w*v4.value.w}, {dx5, value.w*v5.value.w}, {dx6, value.w*v6.value.w}, {dx7, value.w*v7.value.w}, {dx8, value.w*v8.value.w});
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
        __m128 squared_distance = squared_dot_product(this->data.data(), other.data.data(), OutputControl::ALL);
        __m128 distance_sqrt = _mm_sqrt_ps(squared_distance);
        float dist = inv_width*_mm_cvtss_f32(distance_sqrt);
        return EvaluatedResult(dist, value.w*other.value.w);
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

        __m128 distance_sqrt = _mm_sqrt_ps(dist2);
        __m128 distance_binf = _mm_mul_ps(distance_sqrt, _mm_set_ps1(inv_width));       // multiply by the inverse width
        __m128i distance_bin = _mm_cvtps_epi32(distance_binf);                          // cast to int

        __m128 w1 = _mm_set_ps1(value.w);
        __m128 w2 = _mm_set_ps(v4.value.w, v3.value.w, v2.value.w, v1.value.w);
        __m128 weights = _mm_mul_ps(w1, w2);

        QuadEvaluatedResult result;
        _mm_store_si128(reinterpret_cast<__m128i*>(result.distance.data()), distance_bin); // efficient store of bins
        _mm_store_ps(reinterpret_cast<float*>(result.weight.data()), weights);             // efficient store of weights
        return result;
    }

    OctoEvaluatedResult CompactCoordinatesData::evaluate_sse(const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4, const CompactCoordinatesData& v5, const CompactCoordinatesData& v6, const CompactCoordinatesData& v7, const CompactCoordinatesData& v8) const {
        // do two 128-bit calculations
        OctoEvaluatedResult result;

        {
            __m128 dist2_1 = squared_dot_product(this->data.data(), v1.data.data(), OutputControl::FIRST);
            __m128 dist2_2 = squared_dot_product(this->data.data(), v2.data.data(), OutputControl::SECOND);
            __m128 dist2_3 = squared_dot_product(this->data.data(), v3.data.data(), OutputControl::THIRD);
            __m128 dist2_4 = squared_dot_product(this->data.data(), v4.data.data(), OutputControl::FOURTH);

            __m128 dist2 = _mm_add_ps(dist2_1, dist2_2);
            dist2 = _mm_add_ps(dist2, dist2_3);
            dist2 = _mm_add_ps(dist2, dist2_4);

            __m128 distance_sqrt = _mm_sqrt_ps(dist2);
            __m128 distance_binf = _mm_mul_ps(distance_sqrt, _mm_set_ps1(inv_width));
            __m128i distance_bin = _mm_cvtps_epi32(distance_binf);

            __m128 w1 = _mm_set_ps1(value.w);
            __m128 w2 = _mm_set_ps(v4.value.w, v3.value.w, v2.value.w, v1.value.w);
            __m128 weights = _mm_mul_ps(w1, w2);

            _mm_store_si128(reinterpret_cast<__m128i*>(result.distance.data()), distance_bin);
            _mm_store_ps(reinterpret_cast<float*>(result.weight.data()), weights);
        }
        {
            __m128 dist2_1 = squared_dot_product(this->data.data(), v1.data.data(), OutputControl::FIRST);
            __m128 dist2_2 = squared_dot_product(this->data.data(), v2.data.data(), OutputControl::SECOND);
            __m128 dist2_3 = squared_dot_product(this->data.data(), v3.data.data(), OutputControl::THIRD);
            __m128 dist2_4 = squared_dot_product(this->data.data(), v4.data.data(), OutputControl::FOURTH);

            __m128 dist2 = _mm_add_ps(dist2_1, dist2_2);
            dist2 = _mm_add_ps(dist2, dist2_3);
            dist2 = _mm_add_ps(dist2, dist2_4);

            __m128 distance_sqrt = _mm_sqrt_ps(dist2);
            __m128 distance_binf = _mm_mul_ps(distance_sqrt, _mm_set_ps1(inv_width));
            __m128i distance_bin = _mm_cvtps_epi32(distance_binf);

            __m128 w1 = _mm_set_ps1(value.w);
            __m128 w2 = _mm_set_ps(v4.value.w, v3.value.w, v2.value.w, v1.value.w);
            __m128 weights = _mm_mul_ps(w1, w2);

            _mm_store_si128(reinterpret_cast<__m128i*>(result.distance.data()), distance_bin);
            _mm_store_ps(reinterpret_cast<float*>(result.weight.data()) + 4, weights);
        }
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
        __m256 sv12 = _mm256_set_m128(sv1, sv2);    // combine the two 128 bit registers

        __m256 diff = _mm256_sub_ps(svv, sv12);     // calculate the difference

        _mm256_dp_ps(diff, diff, control);          // calculate the dot product, masked to only use the first three components
        return diff;


        // sv[3] = sv1[3] = sv2[3] = 0;                             // zero out the weights
        // __m256 diff = _mm256_sub_ps(svv, sv12);                  // calculate the difference
        // __m256 multiplied = _mm256_mul_ps(diff, diff);           // square the difference

        // __m128 low128 = _mm256_extractf128_ps(multiplied, 0);    // first 128 bits
        // __m128 high128 = _mm256_extractf128_ps(multiplied, 1);   // second 128 bits

        // __m128 sum1 = _mm_hadd_ps(low128, low128);               // sum the first 128 bits
        // __m128 sum2 = _mm_hadd_ps(high128, high128);             // sum the second 128 bits

        // return _mm256_set_m128(sum2, sum1);                      // combine the two sums
    }

    EvaluatedResult CompactCoordinatesData::evaluate_avx(const CompactCoordinatesData& other) const {

    }

    QuadEvaluatedResult CompactCoordinatesData::evaluate_avx(const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4) const {

    }

    OctoEvaluatedResult CompactCoordinatesData::evaluate_avx(const CompactCoordinatesData& v1, const CompactCoordinatesData& v2, const CompactCoordinatesData& v3, const CompactCoordinatesData& v4, const CompactCoordinatesData& v5, const CompactCoordinatesData& v6, const CompactCoordinatesData& v7, const CompactCoordinatesData& v8) const {

    }
#endif