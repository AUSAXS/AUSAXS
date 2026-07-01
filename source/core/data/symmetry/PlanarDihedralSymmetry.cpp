// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <data/symmetry/PlanarDihedralSymmetry.h>

using namespace ausaxs;
using namespace ausaxs::symmetry;

std::unique_ptr<ISymmetry> PlanarDihedralSymmetry::clone() const {return std::make_unique<PlanarDihedralSymmetry>(*this);}

std::span<double> PlanarDihedralSymmetry::span_translation() {
    return std::span<double>(translation.begin(), translation.begin() + 2);
}

Vector3<double> PlanarDihedralSymmetry::group_centre(const Vector3<double>& cm, const Matrix<double>& F) const {
    // (translation.x, translation.y) are in-plane coordinates in the group frame; the axial
    // component is pinned to zero, so F maps the offset into the equatorial plane in world space
    return cm + F*Vector3<double>{translation.x(), translation.y(), 0};
}
