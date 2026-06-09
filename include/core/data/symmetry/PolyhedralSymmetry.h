// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/symmetry/ISymmetry.h>
#include <math/Vector3.h>

namespace ausaxs::symmetry {
    enum class PolyhedralGroup {
        tetrahedral,    //< chiral tetrahedral rotation group T  (12 elements)
        octahedral,     //< chiral octahedral rotation group O   (24 elements)
        icosahedral     //< chiral icosahedral rotation group I  (60 elements)
    };

    /**
     * @brief Chiral polyhedral point-group symmetry (tetrahedral / octahedral / icosahedral).
     *
     * Generates the |G|-1 rotated copies of a body that fill out the full rotation group G.
     * The 12 / 24 / 60 group rotations are fixed; the optimisable parameters are the offset
     * of the body from the group centre and the orientation of the group frame. No mirror
     * operations are included, so all copies remain physically realisable.
     */
    struct PolyhedralSymmetry : public ISymmetry {
        explicit PolyhedralSymmetry(PolyhedralGroup group);

        ISymmetry& add(observer_ptr<const ISymmetry> other) override;
        std::function<Vector3<double>(Vector3<double>)> get_transform(const Vector3<double>& cm, int rep = 1) const override;
        std::unique_ptr<ISymmetry> clone() const override;
        unsigned int repetitions() const override;
        bool is_closed() const override;

        std::span<double> span_translation() override;
        std::span<double> span_rotation() override;
        std::vector<SymmetricDuplicatePair> internal_pair_schedule() const override;

        PolyhedralGroup group;
        Vector3<double> translation{0, 0, 0};   //< offset of the body from the group centre
        Vector3<double> rotation{0, 0, 0};      //< orientation of the group frame (Euler angles)
    };
}
