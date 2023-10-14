#pragma once

#include <math/MathFwd.h>

#include <array>

namespace hist::detail {
    struct EvaluatedResult {
        EvaluatedResult(float distance, float weight) : distance(distance), weight(weight) {}
        float distance, weight;
    };

    struct CompactCoordinatesData {
        CompactCoordinatesData();

        CompactCoordinatesData(const Vector3<double>& v, float w);

        /**
         * @brief Calculate the distance and combined weight between this and another CompactCoordinatesData.
         */
        EvaluatedResult evaluate(const CompactCoordinatesData& other) const;

        union {
            struct {
                float x, y, z, w;
            };
            std::array<float, 4> data;
        };
    };
    static_assert(sizeof(CompactCoordinatesData) == 16, "hist::detail::CompactCoordinates::Data is not 16 bytes");
}