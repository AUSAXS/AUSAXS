#pragma once

#include <math/Vector3.h>

namespace ausaxs::data::detail {
    struct Symmetry {
        Symmetry() = default;
        Symmetry(Vector3<double> translate, Vector3<double> external_axis, double external_angle, Vector3<double> internal_rotate, int repeat) : 
            translate(translate), 
            external_rotate{external_axis, external_angle}, 
            internal_rotate(internal_rotate), 
            repeat(repeat) 
        {}
        Symmetry(Vector3<double> translate) : Symmetry(translate, Vector3<double>(0, 0, 0), 0, Vector3<double>(0, 0, 0), 1) {}
        Symmetry(Vector3<double> external_axis, double angle) : Symmetry(Vector3<double>(0, 0, 0), external_axis, angle, {0, 0, 0}, 1) {}

        // translational vector with respect to the original body
        Vector3<double> translate;

        // external rotation with respect to an axis
        struct {
            Vector3<double> axis;
            double angle;
        } external_rotate;

        // orientation with respect to the original body
        Vector3<double> internal_rotate;

        // the number of times the symmetry should be repeated
        int repeat;
    };
    static_assert(std::is_trivial_v<Symmetry>,          "Symmetry is not trivial");
    static_assert(std::is_standard_layout_v<Symmetry>,  "Symmetry is not standard layout");
    static_assert(supports_nothrow_move_v<Symmetry>,    "Symmetry should support nothrow move semantics.");   
}