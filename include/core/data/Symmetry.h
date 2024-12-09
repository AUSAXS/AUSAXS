#pragma once

#include "math/MatrixUtils.h"
#include <math/Vector3.h>

namespace ausaxs::data::detail {
    struct Symmetry {
        // translational vector with respect to the original body
        Vector3<double> translate = {0, 0, 0};

        // external rotation with respect to an axis
        Vector3<double> external_rotate = {0, 0, 0};
        double external_angle = 0;

        // orientation with respect to the original body
        Vector3<double> internal_rotate = {0, 0, 0};

        // the number of times the symmetry should be repeated
        int repeat = 1;

        template<typename T>
        Vector3<T> apply_transform(const Vector3<T>& v, const Vector3<T>& cm) const {
            auto mi = matrix::rotation_matrix<T>(internal_rotate.x(), internal_rotate.y(), internal_rotate.z());
            auto me = matrix::rotation_matrix<T>(external_rotate, external_angle);
            return me*(mi*(v-cm) + cm + translate);
        }
    };
}