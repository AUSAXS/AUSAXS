// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

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
            assert(initial_relation.translation.dot(repeat_relation.translation) == 0 && "The translation vectors must be orthogonal.");

            // but this is actually not enough, since t_r may also be rotated into the same space as t_i
            // therefore, we must further guarantee that t_r lies in the invariant space of r_i
            assert(
                (
                    repeat_relation.rotation.magnitude() == 0 
                    || matrix::rotation_matrix(repeat_relation.rotation)*repeat_relation.translation == repeat_relation.translation
                )
                && "The translation vector must lie in the invariant space of the rotation matrix."
            );
        }

        Symmetry(_Relation initial_relation, int repetitions = 1)
            : initial_relation(initial_relation), repeat_relation({0, 0, 0}, {0, 0, 0}), repetitions(repetitions) 
        {
            assert(initial_relation.translation.dot(repeat_relation.translation) == 0 && "The translation vectors must be orthogonal.");
            assert(
                (
                    repeat_relation.rotation.magnitude() == 0 
                    || matrix::rotation_matrix(repeat_relation.rotation)*repeat_relation.translation == repeat_relation.translation
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
            _Relation(Vector3<double>&& t, Vector3<double>&& r) : translation(std::move(t)), rotation(std::move(r)) {}
            _Relation(const Vector3<double>& t, const Vector3<double>& r) : translation(t), rotation(r) {}

            Vector3<double> translation;
            Vector3<double> rotation;
            bool operator==(const _Relation& rhs) const = default;
        } initial_relation;

        // The relationship between the N-th repeat and the (N+1)-th repeat
        struct _Repeat {
            _Repeat() = default;
            _Repeat(Vector3<double>&& t, Vector3<double>&& r) : translation(std::move(t)), rotation(std::move(r)) {}
            _Repeat(const Vector3<double>& t, const Vector3<double>& r) : translation(t), rotation(r) {}

            Vector3<double> translation;
            Vector3<double> rotation;
            bool operator==(const _Repeat& rhs) const = default;
        } repeat_relation;
        int repetitions;
    };
    static_assert(std::is_trivial_v<Symmetry>,          "Symmetry is not trivial");
    static_assert(std::is_standard_layout_v<Symmetry>,  "Symmetry is not standard layout");
    static_assert(supports_nothrow_move_v<Symmetry>,    "Symmetry should support nothrow move semantics.");   
}

template<typename Q>
inline std::function<ausaxs::Vector3<Q>(ausaxs::Vector3<Q>)> ausaxs::symmetry::Symmetry::get_transform(const Vector3<double>& cm, int repeat) const {
    Matrix<double>  R_final;
    Vector3<double> T_final;

    {
        auto t_i = initial_relation.translation;
        auto r_i = matrix::rotation_matrix<Q>(initial_relation.rotation);    
        auto t_r = repeat_relation.translation;
        auto r_r = matrix::rotation_matrix<Q>(repeat_relation.rotation);
        assert(t_i.dot(t_r) == 0 && "The translation vectors must be orthogonal.");
        assert(
            (
                repeat_relation.rotation.magnitude() == 0 
                || r_r*t_r == t_r
            )
            && "The translation vector must lie in the invariant space of the rotation matrix."
        );

        // The symmetric ensemble is generated by:
        //   1. Center at origin:                     v  ->  v - cm
        //   2. Place in symmetry space:              v  ->  R_i*(v - cm) + t_i
        //   3. Generate member k by rotating about
        //      the center of symmetry:               s_k = R_r^k * s_0 + sum R_r^j * t_r
        //   4. Transform back to real space so that
        //      member 0 returns to its original 
        //      position (i.e. the original body is
        //      unchanged):                           v_k = R_i^{-1} * (s_k - t_i) + cm
        //
        // Combining:
        //   v_k = R_i^{-1} * R_r^k * R_i * (v - cm) + R_i^{-1} * (R_r^k * t_i - t_i + sum R_r^j * t_r) + cm
        //
        // This can be accumulated iteratively using:
        //   R_conj = R_i^{-1} * R_r * R_i       (conjugation: R_r in the body's frame)
        //   step   = R_i^{-1} * (R_r*t_i + t_r - t_i)
        //   base_T = cm - R_conj*cm + step
        //
        //   R_k = R_conj^k,  T_k = R_conj * T_{k-1} + base_T,  T_1 = base_T

        auto r_i_inv = r_i.T();                                 // R_i^{-1} = R_i^T for rotation matrices
        auto r_conj  = r_i_inv * r_r * r_i;                     // conjugated single-step rotation
        auto step    = r_i_inv * (r_r * t_i + t_r - t_i);       // per-step shift in real space
        auto base_T  = cm - r_conj * cm + step;                 // base translation (k=1)

        R_final = r_conj;
        T_final = base_T;

        // accumulate for repeat > 1
        for (int i = 1; i < repeat; ++i) {
            T_final = r_conj * T_final + base_T;
            R_final = r_conj * R_final;
        }
    }

    if constexpr (std::is_same_v<Q, double>) {
        return [R_final=std::move(R_final), T_final=std::move(T_final)](Vector3<Q> v) {
            return R_final * v + T_final;
        };
    }

    // cast to the correct type for better performance
    Matrix<Q> R_cast = std::move(R_final);
    Vector3<Q> T_cast = std::move(T_final);
    return [R_cast=std::move(R_cast), T_cast=std::move(T_cast)](Vector3<Q> v) {
        return R_cast * v + T_cast;
    };
}

inline bool ausaxs::symmetry::Symmetry::is_closed() const {
    if (repeat_relation.translation.magnitude() != 0) {return false;}

    // due to floating point inaccuracies, we multiply by 100 and round to nearest integer
    // we then compare the value modulo 100*2*pi = 628 with 0
    auto angles = 100*(repetitions+1)*repeat_relation.rotation;
    constexpr int cmp = 100*2*std::numbers::pi;
    return 
        static_cast<int>(std::round(angles.x())) % cmp == 0 && 
        static_cast<int>(std::round(angles.y())) % cmp == 0 && 
        static_cast<int>(std::round(angles.z())) % cmp == 0
    ;
}