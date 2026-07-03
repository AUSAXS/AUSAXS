// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <data/symmetry/PlanarDihedralSymmetry.h>

using namespace ausaxs::symmetry;

template<int N>
std::unique_ptr<ausaxs::symmetry::ISymmetry> PlanarDihedralSymmetry<N>::clone() const {
    return std::make_unique<PlanarDihedralSymmetry<N>>(*this);
}

template<int N>
std::span<double> PlanarDihedralSymmetry<N>::span_translation() {
    return std::span<double>(this->translation.begin(), this->translation.begin() + 2);
}

template<int N>
ausaxs::Vector3<double> PlanarDihedralSymmetry<N>::group_centre(const Vector3<double>& cm, const Matrix<double>& F) const {
    return cm + F*Vector3<double>{this->translation.x(), this->translation.y(), 0};
}

template class ausaxs::symmetry::PlanarDihedralSymmetry<2>;
template class ausaxs::symmetry::PlanarDihedralSymmetry<3>;
template class ausaxs::symmetry::PlanarDihedralSymmetry<4>;
template class ausaxs::symmetry::PlanarDihedralSymmetry<5>;
template class ausaxs::symmetry::PlanarDihedralSymmetry<6>;
template class ausaxs::symmetry::PlanarDihedralSymmetry<7>;
template class ausaxs::symmetry::PlanarDihedralSymmetry<8>;
template class ausaxs::symmetry::PlanarDihedralSymmetry<9>;
template class ausaxs::symmetry::PlanarDihedralSymmetry<10>;
template class ausaxs::symmetry::PlanarDihedralSymmetry<11>;
template class ausaxs::symmetry::PlanarDihedralSymmetry<12>;
