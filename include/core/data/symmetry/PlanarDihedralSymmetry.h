// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/symmetry/DihedralSymmetry.h>

namespace ausaxs::symmetry {
    /**
     * @brief Planar (coplanar) dihedral point group D_n.
     *
     * The same group as DihedralSymmetry but with the body offset constrained to the equatorial plane (perpendicular to the principal axis).
     */
    class PlanarDihedralSymmetry final : public DihedralSymmetry {
        public:
            explicit PlanarDihedralSymmetry(int n) : DihedralSymmetry(n) {}
            std::unique_ptr<ISymmetry> clone() const override;

            // only the two in-plane components are optimisable; the axial component is pinned to zero
            std::span<double> span_translation() override;

        private:
            // interpret the offset in the group frame and drop its axial component, so the copies stay in the equatorial plane regardless of how the frame is oriented
            Vector3<double> group_centre(const Vector3<double>& cm, const Matrix<double>& F) const override;
    };
}
