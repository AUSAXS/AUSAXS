// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/symmetry/DihedralSymmetry.h>

namespace ausaxs::symmetry {
    /**
     * @brief Planar (coplanar) dihedral point group D_n.
     *
     * The same group as DihedralSymmetry — identical elements and pair schedule — but with the body
     * offset constrained to the equatorial plane (perpendicular to the principal axis). A general
     * D_n places the body at some axial height, giving two stacked C_n rings; pinning that height to
     * zero collapses them into a single coplanar ring of 2n copies, the "flat" arrangement that is
     * by far the most common dihedral motif in real assemblies.
     *
     * This removes no pair-distance work (the equal-distance classes depend only on the group
     * rotations, not the offset), but it drops one optimisable degree of freedom, keeping the
     * sampler out of the non-planar states we know a priori are unphysical.
     */
    class PlanarDihedralSymmetry final : public DihedralSymmetry {
        public:
            explicit PlanarDihedralSymmetry(int n) : DihedralSymmetry(n) {}
            std::unique_ptr<ISymmetry> clone() const override;

            // only the two in-plane components are optimisable; the axial component is pinned to zero
            std::span<double> span_translation() override;

        private:
            // interpret the offset in the group frame and drop its axial component, so the copies
            // stay in the equatorial plane regardless of how the frame is oriented
            Vector3<double> group_centre(const Vector3<double>& cm, const Matrix<double>& F) const override;
    };
}
