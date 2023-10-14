#include <hist/detail/CompactCoordinatesData.h>
#include <math/Vector3.h>

#if defined __SSE2__
    #include <nmmintrin.h>
#endif

using namespace hist::detail;

CompactCoordinatesData::CompactCoordinatesData() : x(0), y(0), z(0), w(0) {}

CompactCoordinatesData::CompactCoordinatesData(const Vector3<double>& v, float w) : x(v.x()), y(v.y()), z(v.z()), w(w) {}

constexpr int control_bits = 0b01111111; // only the first three components (x, y, z) are used
EvaluatedResult CompactCoordinatesData::evaluate(const CompactCoordinatesData& other) const {
    #if defined __SSE2__
        // load data into SSE registers
        __m128 this_data = _mm_load_ps(data.data());
        __m128 other_data = _mm_load_ps(other.data.data());
        #if defined __SSE4_1__
            __m128 diff = _mm_sub_ps(this_data, other_data);
            __m128 squared_distance = _mm_dp_ps(diff, diff, control_bits);
        #else
            this_data[3] = other_data[3] = 0; // zero out the fourth component
            __m128 diff = _mm_sub_ps(this_data, other_data);
            __m128 multiplied = _mm_mul_ps(diff, diff);
            __m128 squared_distance = _mm_hadd_ps(multiplied, multiplied);
        #endif
        __m128 distance_sqrt = _mm_sqrt_ps(squared_distance);
        float dist = _mm_cvtss_f32(distance_sqrt);
        return EvaluatedResult(dist, w*other.w);
    #else
        float dx = x - other.x;
        float dy = y - other.y;
        float dz = z - other.z;
        float dist = std::sqrt(dx*dx + dy*dy + dz*dz);
        return EvaluatedResult(dist, w*other.w);
    #endif
}
