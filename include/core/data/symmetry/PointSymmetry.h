// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/symmetry/ISymmetry.h>
#include <math/MatrixUtils.h>

#include <cassert>

namespace ausaxs::symmetry {
    /**
     * @brief A single symmetry operation generating one copy of a body.
     *
     * The copy is produced by applying a fixed translation and rotation relative to the original
     * body. Unlike CyclicSymmetry it is not repeated, so it always contributes exactly one extra copy.
     */
    struct PointSymmetry : public ISymmetry {
        PointSymmetry();
        PointSymmetry(const Vector3<double>& translation, const Vector3<double>& rotation);

        ISymmetry& add(observer_ptr<const ISymmetry> other) override;
        std::unique_ptr<ISymmetry> clone() const override;
        unsigned int repetitions() const override;
        bool is_closed() const override;
        std::string type_name() const override;

        Vector3<double> translation; // Offset of the copy relative to the original body.
        Vector3<double> rotation;    // Rotation of the copy, given as Euler angles.
        std::span<double> span_translation() override;
        std::span<double> span_rotation() override;

    protected:
        AffineTransform _make_transform(const Vector3<double>& anchor, int rep) const override;
    };
}