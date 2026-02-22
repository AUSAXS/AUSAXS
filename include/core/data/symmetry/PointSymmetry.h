// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/symmetry/ISymmetry.h>
#include <math/MatrixUtils.h>

#include <cassert>

namespace ausaxs::symmetry {
    struct PointSymmetry : public ISymmetry {
        PointSymmetry();
        PointSymmetry(const Vector3<double>& translation, const Vector3<double>& rotation);

        ISymmetry& add(observer_ptr<const ISymmetry> other) override;
        std::function<Vector3<double>(Vector3<double>)> get_transform(const Vector3<double>& cm, int rep = 1) const override;
        std::unique_ptr<ISymmetry> clone() const override;
        unsigned int repetitions() const override;
        bool is_closed() const override;

        Vector3<double> translation;
        Vector3<double> rotation;
        std::span<double> span_translation() override;
        std::span<double> span_rotation() override;
    };
}