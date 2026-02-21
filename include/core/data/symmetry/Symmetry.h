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
    struct Symmetry : public ISymmetry {
        Symmetry();
        Symmetry(_Relation initial_relation, _Repeat repeat_relation, int repetitions = 1);
        Symmetry(Vector3<double> offset, Vector3<double> repeat_translation, Vector3<double> repeat_axis, double repeat_rotation, int repetitions = 1);

        ISymmetry& add(observer_ptr<const ISymmetry> other) override;
        std::function<Vector3<double>(Vector3<double>)> get_transform(const Vector3<double>& cm, int rep = 1) const override;
        std::unique_ptr<ISymmetry> clone() const override;
        bool is_closed() const override;
    };
}