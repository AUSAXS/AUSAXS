// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <data/symmetry/DihedralSymmetry.h>
#include <math/MatrixUtils.h>

#include <numbers>

using namespace ausaxs::symmetry;

template<int N>
std::unique_ptr<ausaxs::symmetry::ISymmetry> DihedralSymmetry<N>::clone() const {
    return std::make_unique<DihedralSymmetry<N>>(*this);
}

template<int N>
const IPolyhedralSymmetry::GroupData& DihedralSymmetry<N>::group() const {
    static const GroupData data = build({
        matrix::rotation_matrix<double>({0, 0, 1}, 2*std::numbers::pi/N),
        matrix::rotation_matrix<double>({1, 0, 0}, std::numbers::pi)
    }, 2*N);
    return data;
}

// only these orders are reachable through the symmetry::type enum / name parser
template class ausaxs::symmetry::DihedralSymmetry<2>;
template class ausaxs::symmetry::DihedralSymmetry<3>;
template class ausaxs::symmetry::DihedralSymmetry<4>;
template class ausaxs::symmetry::DihedralSymmetry<5>;
template class ausaxs::symmetry::DihedralSymmetry<6>;
template class ausaxs::symmetry::DihedralSymmetry<7>;
template class ausaxs::symmetry::DihedralSymmetry<8>;
template class ausaxs::symmetry::DihedralSymmetry<9>;
template class ausaxs::symmetry::DihedralSymmetry<10>;
template class ausaxs::symmetry::DihedralSymmetry<11>;
template class ausaxs::symmetry::DihedralSymmetry<12>;
