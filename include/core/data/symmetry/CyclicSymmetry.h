// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/symmetry/ISymmetry.h>
#include <math/MatrixUtils.h>

#include <cassert>

namespace ausaxs::symmetry {
    /**
     * @brief Standard rotational symmetry (c2, c3, ...) parameterised as centre-of-rotation offset + axis/angle.
     *
     * The optimisable parameters are:
     *   initial_relation.translation — offset of the original body from the rotation centre.
     *   repeat_relation.axis / angle — per-step rotation for generating each copy.
     *   repeat_relation.translation  — optional per-step axial translation (screw symmetries).
     */
    struct CyclicSymmetry : public ISymmetry {
        struct _Relation {
            _Relation() : translation{0, 0, 0} {}
            _Relation(const Vector3<double>& t) : translation(t) {}
            Vector3<double> translation;
            bool operator==(const _Relation&) const = default;
        };

        struct _Repeat {
            _Repeat() : translation{0, 0, 0}, axis{0, 0, 1}, angle(0) {}
            _Repeat(const Vector3<double>& t, const Vector3<double>& ax, double ang)
                : translation(t), axis(ax), angle(ang) {}
            _Repeat(const Vector3<double>& ax, double ang)
                : translation{0, 0, 0}, axis(ax), angle(ang) {}
            Vector3<double> translation;
            Vector3<double> axis;
            double          angle = 0;
            bool operator==(const _Repeat&) const = default;
        };

        _Relation _initial_relation;
        _Repeat   _repeat_relation;
        int       _repetitions = 1;

        CyclicSymmetry();
        CyclicSymmetry(_Relation initial_relation, _Repeat repeat_relation, int repetitions = 1);
        CyclicSymmetry(Vector3<double> offset, Vector3<double> repeat_translation, Vector3<double> repeat_axis, double repeat_rotation, int repetitions = 1);

        ISymmetry& add(observer_ptr<const ISymmetry> other) override;
        std::function<Vector3<double>(Vector3<double>)> get_transform(const Vector3<double>& cm, int rep = 1) const override;
        std::unique_ptr<ISymmetry> clone() const override;
        unsigned int repetitions() const override;
        bool is_closed() const override;

        std::span<double> span_translation() override;
        std::span<double> span_rotation() override;
    };
}