#pragma once

#include <math/Vector3.h>
#include <math/MatrixUtils.h>

namespace ausaxs::data::detail {
    struct Symmetry {
        Symmetry() = default;
        Symmetry(Vector3<double> translate, Vector3<double> external_axis, double external_angle, Vector3<double> internal_rotate, int repeat) : 
            translate(translate), 
            external_rotate{external_axis, external_angle}, 
            internal_rotate(internal_rotate), 
            repeat(repeat) 
        {}
        Symmetry(Vector3<double> translate, Vector3<double> external_axis, double external_angle, Vector3<double> internal_rotate) : 
            Symmetry(translate, external_axis, external_angle, internal_rotate, 1) 
        {}
        Symmetry(Vector3<double> translate) : Symmetry(translate, Vector3<double>(0, 0, 0), 0, Vector3<double>(0, 0, 0), 1) {}
        Symmetry(Vector3<double> external_axis, double angle) : Symmetry(Vector3<double>(0, 0, 0), external_axis, angle, {0, 0, 0}, 1) {}

        /**
         * @brief Get the transformation function for this symmetry.
         */
        template<typename T>
        std::function<Vector3<T>(Vector3<T>)> get_transform(Vector3<T> cm) const;

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

template<typename T>
std::function<ausaxs::Vector3<T>(ausaxs::Vector3<T>)> ausaxs::data::detail::Symmetry::get_transform(Vector3<T> cm) const {
    bool translate = this->translate != Vector3<T>(0, 0, 0);
    bool internal_rotate = this->internal_rotate != Vector3<T>(0, 0, 0);
    bool external_rotate = this->external_rotate.angle != 0;

    if (internal_rotate && external_rotate) {
        Matrix<T> internal_rotate = matrix::rotation_matrix<T>(this->internal_rotate);
        Matrix<T> external_rotate = matrix::rotation_matrix<T>(this->external_rotate.axis, this->external_rotate.angle);

        if (translate) {
            return [cm, r_int=internal_rotate, r_ext=external_rotate, t=this->translate](Vector3<T> v) {
                return r_ext*(r_int*(v-cm) + cm) + t;
            };
        } else {
            return [cm, r_int=internal_rotate, r_ext=external_rotate](Vector3<T> v) {
                return r_ext*(r_int*(v-cm) + cm);
            };
        }
    }

    if (internal_rotate) {
        Matrix<T> internal_rotate = matrix::rotation_matrix<T>(this->internal_rotate);

        if (translate) {
            return [cm, r_int=internal_rotate, t=this->translate](Vector3<T> v) {
                return r_int*(v-cm) + cm + t;
            };
        } else {
            return [cm, r_int=internal_rotate](Vector3<T> v) {
                return r_int*(v-cm) + cm;
            };
        }
    }

    if (external_rotate) {
        Matrix<T> external_rotate = matrix::rotation_matrix<T>(this->external_rotate.axis, this->external_rotate.angle);

        if (translate) {
            return [r_ext=external_rotate, t=this->translate](Vector3<T> v) {
                return r_ext*v + t;
            };
        } else {
            return [r_ext=external_rotate](Vector3<T> v) {
                return r_ext*v;
            };
        }
    }

    if (translate) {
        return [t=this->translate](Vector3<T> v) {
            return v + t;
        };
    }

    return [] (Vector3<T> v) {return v;};
}