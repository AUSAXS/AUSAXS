#pragma once

#include <math/Vector3.h>
#include <math/MatrixUtils.h>

#include <numbers>
#include <cassert>

namespace ausaxs::symmetry {
    struct Symmetry {
        struct _Relation; struct _Repeat;

        Symmetry() = default;
        Symmetry(_Relation initial_relation, _Repeat repeat_relation, int repetitions = 1) 
            : initial_relation(initial_relation), repeat_relation(repeat_relation), repetitions(repetitions)
        {
            // we must guarantee that all symmetric duplicates are equidistant, which is complicated because we have two separate translation vectors
            // if t_i and t_r are linearly dependent, then the distance to the first repeat will be offset by t_i compared to the following duplicates
            // therefore, the action of t_i and t_r cannot overlap, i.e. they must be orthogonal
            assert(initial_relation.translation.dot(repeat_relation.translate) == 0 && "The translation vectors must be orthogonal.");

            // but this is actually not enough, since t_r may also be rotated into the same space as t_i
            // therefore, we must further guarantee that t_r lies in the invariant space of r_i
            assert(
                (
                    repeat_relation.rotate.magnitude() == 0 
                    || matrix::rotation_matrix(repeat_relation.rotate)*repeat_relation.translate == repeat_relation.translate
                )
                && "The translation vector must lie in the invariant space of the rotation matrix."
            );
        }

        Symmetry(_Relation initial_relation, int repetitions = 1)
            : initial_relation(initial_relation), repeat_relation({0, 0, 0}, {0, 0, 0}), repetitions(repetitions) 
        {
            assert(initial_relation.translation.dot(repeat_relation.translate) == 0 && "The translation vectors must be orthogonal.");
            assert(
                (
                    repeat_relation.rotate.magnitude() == 0 
                    || matrix::rotation_matrix(repeat_relation.rotate)*repeat_relation.translate == repeat_relation.translate
                )
                && "The translation vector must lie in the invariant space of the rotation matrix."
            );
        }

        /**
         * @brief Get the transform taking the original structure to the N-th repeat. 
         *
         * @param cm The center of mass of the original structure.
         * @param repeat The number of times the symmetry should be repeated.
         */
        template<typename T>
        std::function<Vector3<T>(Vector3<T>)> get_transform(const Vector3<double>& cm, int repeat = 1) const;

        /**
         * @brief Determine if the symmetry is closed, i.e. the repeat+1-th transformation is the identity.
         */
        bool is_closed() const;

        bool operator==(const Symmetry& rhs) const = default;

        // The relationship between the original body and the first repeat of this symmetry
        struct _Relation {
            _Relation() = default;
            _Relation(Vector3<double>&& t, Vector3<double>&& o) : translation(std::move(t)), orientation(std::move(o)) {}
            _Relation(const Vector3<double>& t, const Vector3<double>& o) : translation(t), orientation(o) {}

            Vector3<double> translation;
            Vector3<double> orientation;
            bool operator==(const _Relation& rhs) const = default;
        } initial_relation;

        // The relationship between the N-th repeat and the (N+1)-th repeat
        struct _Repeat {
            _Repeat() = default;
            _Repeat(Vector3<double>&& t, Vector3<double>&& r) : translate(std::move(t)), rotate(std::move(r)) {}
            _Repeat(const Vector3<double>& t, const Vector3<double>& r) : translate(t), rotate(r) {}

            Vector3<double> translate;
            Vector3<double> rotate;
            bool operator==(const _Repeat& rhs) const = default;
        } repeat_relation;
        int repetitions;
    };
    static_assert(std::is_trivial_v<Symmetry>,          "Symmetry is not trivial");
    static_assert(std::is_standard_layout_v<Symmetry>,  "Symmetry is not standard layout");
    static_assert(supports_nothrow_move_v<Symmetry>,    "Symmetry should support nothrow move semantics.");   
}

template<typename Q>
inline std::function<ausaxs::Vector3<Q>(ausaxs::Vector3<Q>)> ausaxs::symmetry::Symmetry::get_transform(const Vector3<double>& t_cm, int repeat) const {
    Matrix<double>  R_final;
    Vector3<double> T_final;

    {
        auto t_i = initial_relation.translation;
        auto r_i = matrix::rotation_matrix<Q>(initial_relation.orientation);    
        auto t_r = repeat_relation.translate;
        auto r_r = matrix::rotation_matrix<Q>(repeat_relation.rotate);
        assert(t_i.dot(t_r) == 0 && "The translation vectors must be orthogonal.");
        assert(
            (
                repeat_relation.rotate.magnitude() == 0 
                || r_r*t_r == t_r
            )
            && "The translation vector must lie in the invariant space of the rotation matrix."
        );

        // initial transform
        // v' = R_r * (R_i * (v - t_cm) + t_cm + t_i) + t_r
        //    = R_r*R_i*v - R_r*R_i*t_cm + R_r*t_cm + R_r*t_i + t_r
        //    = R_r*R_i*v + (R_r*t_cm - R_r*R_i*t_cm + R_r*t_i + t_r)
        //
        // subsequent transforms
        // v' = R_r * v + t_r

        R_final = r_r * r_i;
        T_final = r_r*(t_cm + t_i - r_i*t_cm) + t_r;

        // accumulate transformations from 1 to repeat
        for (int i = 1; i < repeat; ++i) {
            R_final = r_r*R_final;          // accumulate the rotation
            T_final = r_r*T_final + t_r;    // accumulate the translation
        }
    }

    if constexpr (std::is_same_v<Q, double>) {
        return [R_final=std::move(R_final), T_final=std::move(T_final)](Vector3<Q> v) {
            return R_final * v + T_final;
        };
    }

    // cast to the correct type for better performance
    Matrix<Q> R_cast = R_final;
    Vector3<Q> T_cast = T_final;
    return [R_cast=std::move(R_cast), T_cast=std::move(T_cast)](Vector3<Q> v) {
        return R_cast * v + T_cast;
    };
}

inline bool ausaxs::symmetry::Symmetry::is_closed() const {
    if (repeat_relation.translate.magnitude() != 0) {return false;}

    // due to floating point inaccuracies, we multiply by 100 and round to nearest integer
    // we then compare the value modulo 100*2*pi = 628 with 0
    auto angles = 100*(repetitions+1)*repeat_relation.rotate;
    constexpr int cmp = 100*2*std::numbers::pi;
    return 
        static_cast<int>(std::round(angles.x())) % cmp == 0 && 
        static_cast<int>(std::round(angles.y())) % cmp == 0 && 
        static_cast<int>(std::round(angles.z())) % cmp == 0
    ;
}