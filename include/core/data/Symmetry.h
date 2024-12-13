#pragma once

#include <math/Vector3.h>
#include <math/MatrixUtils.h>

namespace ausaxs::data::detail {
    struct Symmetry {
        Symmetry() = default;
        Symmetry(Vector3<double> translate, Vector3<double> center, Vector3<double> angles, Vector3<double> internal_rotate, int repeat) : 
            translate(translate), 
            external_rotate{center, angles}, 
            internal_rotate(internal_rotate), 
            repeat(repeat) 
        {}
        Symmetry(Vector3<double> translate, Vector3<double> center, Vector3<double> angles, Vector3<double> internal_rotate) : 
            Symmetry(translate, center, angles, internal_rotate, 1) 
        {}
        Symmetry(Vector3<double> translate) : Symmetry(translate, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, 1) {}
        Symmetry(Vector3<double> center, Vector3<double> angles) : Symmetry({0, 0, 0}, center, angles, {0, 0, 0}, 1) {}

        /**
         * @brief Get the transformation function for this symmetry.
         */
        template<typename T>
        std::function<Vector3<T>(Vector3<T>)> get_transform(Vector3<T> cm) const;

        // translational vector with respect to the original body
        Vector3<double> translate;

        // external rotation with respect to an axis
        struct {
            Vector3<double> center;
            Vector3<double> angle;
        } external_rotate;

        // orientation with respect to the original body
        struct {
            Vector3<double> angle;
        } internal_rotate;

        // the number of times the symmetry should be repeated
        int repeat;
    };
    static_assert(std::is_trivial_v<Symmetry>,          "Symmetry is not trivial");
    static_assert(std::is_standard_layout_v<Symmetry>,  "Symmetry is not standard layout");
    static_assert(supports_nothrow_move_v<Symmetry>,    "Symmetry should support nothrow move semantics.");   
}

template<typename Q>
std::function<ausaxs::Vector3<Q>(ausaxs::Vector3<Q>)> ausaxs::data::detail::Symmetry::get_transform(Vector3<Q> cm) const {
    // v' = r_ext * (r_int * (v - p_cm) + p_cm - p_sym) + p_sym + t
    //    = r_ext*r_int*v - r_ext*r_int*p_cm + r_ext*p_cm - r_ext*p_sym + p_sym + t
    //    = r_ext*r_int*v + (r_ext*p_cm - r_ext*r_int*p_cm - r_ext*p_sym + p_sym + t)
    auto r_ext = matrix::rotation_matrix<Q>(external_rotate.angle);
    auto r_int = matrix::rotation_matrix<Q>(internal_rotate.angle);
    auto p_sym = external_rotate.center;
    auto p_cm = cm;
    auto t = translate;

    Matrix<Q> R = r_ext*r_int;
    Vector3<Q> T = r_ext*p_cm - R*p_cm- r_ext*p_sym + p_sym + t;
    return [R, T] (Vector3<Q> v) {return R*v + T;};
}