// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/symmetry/ISymmetry.h>
#include <math/MatrixUtils.h>

#include <cassert>

namespace ausaxs::symmetry {
    /**
     * @brief Decoupled dimer-type symmetry where translation and orientation are independent parameters.
     *
     * In contrast to the standard Symmetry class — where the offset and orientation are geometrically
     * linked through a shared rotation centre — PointSymmetry exposes them independently:
     *
     *   initial_relation.translation  = d   : direct CM-to-CM offset vector (lab frame).
     *   repeat_relation.axis / angle       : orientation of the copy (R = rotation_matrix(axis, angle)).
     *
     * The transform for the single copy (rep must be 1) is:
     *   v' = R(axis, angle) * (v - cm) + cm + d
     *      = R * v + (I - R)*cm + d
     *
     * This parameterisation fully decouples translation from orientation, making it much easier
     * for an optimiser to converge for dimer structures where both need to be searched independently.
     *
     * repeat_relation.translation is unused; repetitions is always 1.
     */
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