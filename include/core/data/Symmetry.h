#pragma once

#include <math/Vector3.h>
#include <math/MatrixUtils.h>

namespace ausaxs::data::detail {
    struct Symmetry {
        Symmetry() = default;
        Symmetry(Vector3<double> translate, Vector3<double> center, Vector3<double> angles, int repeat = 1) : 
            translate(translate), 
            external_rotate{center, angles}, 
            repeat(repeat) 
        {}
        Symmetry(Vector3<double> translate) : Symmetry(translate, {0, 0, 0}, {0, 0, 0}, 1) {}
        Symmetry(Vector3<double> center, Vector3<double> angles) : Symmetry({0, 0, 0}, center, angles, 1) {}

        /**
         * @brief Get the transform taking the original structure to the N-th repeat. 
         *
         * @param cm The center of mass of the original structure.
         * @param repeat The number of times the symmetry should be repeated.
         */
        template<typename T>
        std::function<Vector3<T>(Vector3<T>)> get_transform(int repeat = 1) const;

        bool operator==(const Symmetry& rhs) const = default;

        // translational vector with respect to the original body
        Vector3<double> translate;

        // external rotation with respect to an axis
        struct _extrot {
            Vector3<double> center;
            Vector3<double> angle;
            bool operator==(const _extrot& rhs) const = default;
        } external_rotate;

        // the number of times the symmetry should be repeated
        int repeat;
    };
    static_assert(std::is_trivial_v<Symmetry>,          "Symmetry is not trivial");
    static_assert(std::is_standard_layout_v<Symmetry>,  "Symmetry is not standard layout");
    static_assert(supports_nothrow_move_v<Symmetry>,    "Symmetry should support nothrow move semantics.");   
}

template<typename Q>
std::function<ausaxs::Vector3<Q>(ausaxs::Vector3<Q>)> ausaxs::data::detail::Symmetry::get_transform(int repeat) const {
    // accumulate transformations from 1 to repeat
    Matrix<double>  R_final = matrix::identity(3);
    Vector3<double> T_final(0, 0, 0);
    for (int i = 0; i < repeat; ++i) {
        // transformation incorporating translation, rotation about an axis, and orientation (outdated)
        // v' = r_ext * (r_int * (v - p_cm) + p_cm - p_sym) + p_sym + t
        //    = r_ext*r_int*v - r_ext*r_int*p_cm + r_ext*p_cm - r_ext*p_sym + p_sym + t
        //    = r_ext*r_int*v + (r_ext*p_cm - r_ext*r_int*p_cm - r_ext*p_sym + p_sym + t)

        // transformation incorporating translation and rotation about an axis
        // v' = r_ext * (v - p_sym) + p_sym + t
        //    = r_ext*v - r_ext*p_sym + p_sym + t
        auto r_ext = matrix::rotation_matrix<Q>(external_rotate.angle);
        auto p_sym = external_rotate.center;
        auto t = translate;

        Matrix<Q>  R = r_ext;
        Vector3<Q> T = r_ext*p_sym + p_sym + t;

        R_final = R*R_final;        // accumulate the rotation
        T_final = R*T_final + T;    // accumulate the translation
    }

    if constexpr (std::is_same_v<Q, double>) {
        return [R_final, T_final](Vector3<Q> v) {
            return R_final * v + T_final;
        };
    } 

    // cast to the correct type for better performance
    Matrix<Q> R_cast = R_final;
    Vector3<Q> T_cast = T_final;
    return [R_cast, T_cast](Vector3<Q> v) {
        return R_cast * v + T_cast;
    };    
}