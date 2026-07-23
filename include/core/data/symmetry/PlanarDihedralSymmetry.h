// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/symmetry/DihedralSymmetry.h>

namespace ausaxs::symmetry {
    /**
     * @brief Planar (coplanar) dihedral point group D_N.
     */
    template<int N>
    class PlanarDihedralSymmetry final : public DihedralSymmetry<N> {
        public:
            std::unique_ptr<ISymmetry> clone() const override;
            std::string type_name() const override;

            std::span<double> span_translation() override;

        private:
            Vector3<double> group_centre(const Vector3<double>& cm, const Matrix<double>& F) const override;
    };
}
